//! Core data structures for the dual graph system.
//! 
//! This module implements the Arena pattern with index-based references
//! to avoid Rc<RefCell<T>> and handle circular references cleanly.
//!
//! The dual graph consists of:
//! - **Center** (Region): Voronoi polygon, corresponds to an original input point
//! - **Corner**: Voronoi vertex, corresponds to a Delaunay triangle circumcenter
//! - **Edge**: Connects two Centers (Delaunay edge) and two Corners (Voronoi edge)

use crate::geometry::Point;

/// Sentinel value for "no reference" (like null).
pub const NONE: usize = usize::MAX;

/// Biome types determined by elevation and moisture.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Biome {
    #[default]
    Ocean,
    Lake,
    Ice,
    Marsh,
    Beach,
    Snow,
    Tundra,
    Bare,
    Scorched,
    Taiga,
    Shrubland,
    TemperateDesert,
    TemperateRainForest,
    TemperateDeciduousForest,
    Grassland,
    TropicalRainForest,
    TropicalSeasonalForest,
    SubtropicalDesert,
}

/// A Voronoi polygon center (Region).
/// Corresponds to one of the original input points.
#[derive(Debug, Clone, Default)]
pub struct Center {
    /// Index in the centers array.
    pub index: usize,
    /// Position (the original input point, or centroid after Lloyd relaxation).
    pub position: Point,

    // === Topology (indices) ===
    /// Neighboring Center indices (share an edge).
    pub neighbors: Vec<usize>,
    /// Border Edge indices.
    pub borders: Vec<usize>,
    /// Corner indices (vertices of the Voronoi polygon).
    pub corners: Vec<usize>,

    // === Terrain properties ===
    /// True if this is water (ocean or lake).
    pub water: bool,
    /// True if this is ocean (water connected to map boundary).
    pub ocean: bool,
    /// True if this is a coast (land adjacent to ocean).
    pub coast: bool,
    /// Biome type.
    pub biome: Biome,
    /// Elevation (0.0 = sea level, 1.0 = mountain peak).
    pub elevation: f64,
    /// Moisture (0.0 = dry, 1.0 = wet).
    pub moisture: f64,
}

impl Center {
    pub fn new(index: usize, position: Point) -> Self {
        Self {
            index,
            position,
            ..Default::default()
        }
    }
}

/// A Voronoi vertex (Corner).
/// Corresponds to the circumcenter of a Delaunay triangle.
#[derive(Debug, Clone, Default)]
pub struct Corner {
    /// Index in the corners array.
    pub index: usize,
    /// Position (circumcenter of the Delaunay triangle).
    pub position: Point,

    // === Topology (indices) ===
    /// Adjacent Corner indices (connected by Voronoi edges).
    pub adjacent: Vec<usize>,
    /// Edge indices that protrude from this corner.
    pub protrudes: Vec<usize>,
    /// Center indices that this corner touches (vertices of).
    pub touches: Vec<usize>,

    // === Terrain properties ===
    /// True if this corner is water.
    pub water: bool,
    /// True if this corner is in the ocean.
    pub ocean: bool,
    /// True if this is a coast corner.
    pub coast: bool,
    /// Elevation at this corner.
    pub elevation: f64,
    /// Moisture at this corner.
    pub moisture: f64,
    /// River volume flowing through this corner.
    pub river: i32,
    /// Index of the downslope corner (for river flow), NONE if no downslope.
    pub downslope: usize,
    /// Index of the watershed corner (coast corner this flows to).
    pub watershed: usize,
    /// Size of the watershed.
    pub watershed_size: i32,

    // === Boundary flags ===
    /// True if this corner is on the map boundary.
    pub border: bool,
}

impl Corner {
    pub fn new(index: usize, position: Point) -> Self {
        Self {
            index,
            position,
            downslope: NONE,
            watershed: NONE,
            ..Default::default()
        }
    }
}

/// An edge in the dual graph.
/// Connects two Centers (Delaunay edge) and two Corners (Voronoi edge).
#[derive(Debug, Clone, Default)]
pub struct Edge {
    /// Index in the edges array.
    pub index: usize,

    // === Delaunay edge (connects two Centers) ===
    /// Center on one side, NONE if ghost/boundary.
    pub d0: usize,
    /// Center on the other side, NONE if ghost/boundary.
    pub d1: usize,

    // === Voronoi edge (connects two Corners) ===
    /// Corner at one end, NONE if extends to infinity.
    pub v0: usize,
    /// Corner at the other end, NONE if extends to infinity.
    pub v1: usize,

    // === Properties ===
    /// Midpoint of the edge.
    pub midpoint: Point,
    /// River volume flowing along this edge.
    pub river: i32,
}

impl Edge {
    pub fn new(index: usize) -> Self {
        Self {
            index,
            d0: NONE,
            d1: NONE,
            v0: NONE,
            v1: NONE,
            ..Default::default()
        }
    }
}

/// The complete dual mesh structure holding all graph elements.
/// Uses Arena pattern - all references are usize indices into these vectors.
#[derive(Debug, Clone, Default)]
pub struct DualMesh {
    /// All Centers (Voronoi polygons / regions).
    pub centers: Vec<Center>,
    /// All Corners (Voronoi vertices).
    pub corners: Vec<Corner>,
    /// All Edges.
    pub edges: Vec<Edge>,

    // === Ghost element tracking ===
    /// Number of "solid" (non-ghost) centers.
    pub num_solid_centers: usize,
    /// Number of "solid" (non-ghost) corners.
    pub num_solid_corners: usize,
    /// Number of "solid" (non-ghost) edges.
    pub num_solid_edges: usize,
}

impl DualMesh {
    /// Create a new empty dual mesh.
    pub fn new() -> Self {
        Self::default()
    }

    // === Ghost element checks ===

    /// Check if a center index refers to a ghost element.
    #[inline]
    pub fn is_ghost_center(&self, index: usize) -> bool {
        index >= self.num_solid_centers
    }

    /// Check if a corner index refers to a ghost element.
    #[inline]
    pub fn is_ghost_corner(&self, index: usize) -> bool {
        index >= self.num_solid_corners
    }

    /// Check if an edge index refers to a ghost element.
    #[inline]
    pub fn is_ghost_edge(&self, index: usize) -> bool {
        index >= self.num_solid_edges
    }

    // === Accessors ===

    /// Get a center by index.
    #[inline]
    pub fn center(&self, index: usize) -> Option<&Center> {
        self.centers.get(index)
    }

    /// Get a mutable center by index.
    #[inline]
    pub fn center_mut(&mut self, index: usize) -> Option<&mut Center> {
        self.centers.get_mut(index)
    }

    /// Get a corner by index.
    #[inline]
    pub fn corner(&self, index: usize) -> Option<&Corner> {
        self.corners.get(index)
    }

    /// Get a mutable corner by index.
    #[inline]
    pub fn corner_mut(&mut self, index: usize) -> Option<&mut Corner> {
        self.corners.get_mut(index)
    }

    /// Get an edge by index.
    #[inline]
    pub fn edge(&self, index: usize) -> Option<&Edge> {
        self.edges.get(index)
    }

    /// Get a mutable edge by index.
    #[inline]
    pub fn edge_mut(&mut self, index: usize) -> Option<&mut Edge> {
        self.edges.get_mut(index)
    }

    // === Iteration helpers ===

    /// Iterate over all solid (non-ghost) centers.
    pub fn solid_centers(&self) -> impl Iterator<Item = &Center> {
        self.centers.iter().take(self.num_solid_centers)
    }

    /// Iterate over all solid (non-ghost) corners.
    pub fn solid_corners(&self) -> impl Iterator<Item = &Corner> {
        self.corners.iter().take(self.num_solid_corners)
    }

    /// Iterate over all solid (non-ghost) edges.
    pub fn solid_edges(&self) -> impl Iterator<Item = &Edge> {
        self.edges.iter().take(self.num_solid_edges)
    }

    /// Get the Voronoi polygon vertices for a center (in order).
    /// Optionally clamps vertices to the given bounds.
    pub fn get_polygon_vertices_clipped(&self, center_index: usize, bounds: Option<&crate::geometry::BoundingBox>) -> Vec<Point> {
        let center = match self.center(center_index) {
            Some(c) => c,
            None => return vec![],
        };

        // Collect corner positions, optionally clamping to bounds
        let mut vertices: Vec<Point> = center.corners
            .iter()
            .filter_map(|&ci| {
                self.corner(ci).map(|c| {
                    if let Some(b) = bounds {
                        // Clamp to bounds
                        Point::new(
                            c.position.x.clamp(b.min.x, b.max.x),
                            c.position.y.clamp(b.min.y, b.max.y),
                        )
                    } else {
                        c.position
                    }
                })
            })
            .collect();

        if vertices.len() < 3 {
            return vertices;
        }

        // Sort corners by angle around center
        let cx = center.position.x;
        let cy = center.position.y;
        vertices.sort_by(|a, b| {
            let angle_a = (a.y - cy).atan2(a.x - cx);
            let angle_b = (b.y - cy).atan2(b.x - cx);
            angle_a.partial_cmp(&angle_b).unwrap_or(std::cmp::Ordering::Equal)
        });

        vertices
    }

    /// Get the Voronoi polygon vertices for a center (in order).
    pub fn get_polygon_vertices(&self, center_index: usize) -> Vec<Point> {
        self.get_polygon_vertices_clipped(center_index, None)
    }

    // === Debug validation ===

    /// Validate the mesh structure (debug builds only).
    #[cfg(debug_assertions)]
    pub fn validate(&self) -> Result<(), String> {
        // Check center references
        for center in &self.centers {
            for &neighbor in &center.neighbors {
                if neighbor != NONE && neighbor >= self.centers.len() {
                    return Err(format!(
                        "Center {} has invalid neighbor index {}",
                        center.index, neighbor
                    ));
                }
            }
            for &corner in &center.corners {
                if corner != NONE && corner >= self.corners.len() {
                    return Err(format!(
                        "Center {} has invalid corner index {}",
                        center.index, corner
                    ));
                }
            }
            for &border in &center.borders {
                if border != NONE && border >= self.edges.len() {
                    return Err(format!(
                        "Center {} has invalid border index {}",
                        center.index, border
                    ));
                }
            }
        }

        // Check corner references
        for corner in &self.corners {
            for &adjacent in &corner.adjacent {
                if adjacent != NONE && adjacent >= self.corners.len() {
                    return Err(format!(
                        "Corner {} has invalid adjacent index {}",
                        corner.index, adjacent
                    ));
                }
            }
            for &touches in &corner.touches {
                if touches != NONE && touches >= self.centers.len() {
                    return Err(format!(
                        "Corner {} has invalid touches index {}",
                        corner.index, touches
                    ));
                }
            }
        }

        // Check edge references
        for edge in &self.edges {
            if edge.d0 != NONE && edge.d0 >= self.centers.len() {
                return Err(format!("Edge {} has invalid d0 index {}", edge.index, edge.d0));
            }
            if edge.d1 != NONE && edge.d1 >= self.centers.len() {
                return Err(format!("Edge {} has invalid d1 index {}", edge.index, edge.d1));
            }
            if edge.v0 != NONE && edge.v0 >= self.corners.len() {
                return Err(format!("Edge {} has invalid v0 index {}", edge.index, edge.v0));
            }
            if edge.v1 != NONE && edge.v1 >= self.corners.len() {
                return Err(format!("Edge {} has invalid v1 index {}", edge.index, edge.v1));
            }
        }

        Ok(())
    }

    #[cfg(not(debug_assertions))]
    pub fn validate(&self) -> Result<(), String> {
        Ok(())
    }
}
