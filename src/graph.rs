//! Build the dual graph (Voronoi/Delaunay) from delaunator output.
//!
//! This is the critical step that transforms raw triangulation data
//! into our Arena-based DualMesh structure.

use std::collections::{HashMap, HashSet};

use crate::geometry::{Point, BoundingBox, circumcenter};
use crate::mesh::{DualMesh, Center, Corner, Edge, NONE};

/// Build a DualMesh from a set of points.
pub fn build_dual_mesh(points: &[Point], bounds: &BoundingBox) -> DualMesh {
    if points.len() < 3 {
        return DualMesh::new();
    }

    // Convert to delaunator points
    let delaunator_points: Vec<delaunator::Point> = points
        .iter()
        .map(|p| delaunator::Point { x: p.x, y: p.y })
        .collect();

    // Perform Delaunay triangulation
    let triangulation = delaunator::triangulate(&delaunator_points);
    if triangulation.triangles.is_empty() {
        return DualMesh::new();
    }

    let mut mesh = DualMesh::new();

    // Step 1: Create Centers (one per input point)
    for (i, point) in points.iter().enumerate() {
        let mut center = Center::new(i, *point);
        // Mark boundary points
        if bounds.is_on_boundary(point, 1.0) {
            center.water = true;
            center.ocean = true;
        }
        mesh.centers.push(center);
    }
    mesh.num_solid_centers = mesh.centers.len();

    // Step 2: Create Corners (one per triangle, at circumcenter)
    let num_triangles = triangulation.triangles.len() / 3;
    for t in 0..num_triangles {
        let i0 = triangulation.triangles[3 * t];
        let i1 = triangulation.triangles[3 * t + 1];
        let i2 = triangulation.triangles[3 * t + 2];

        let p0 = &points[i0];
        let p1 = &points[i1];
        let p2 = &points[i2];

        let cc = circumcenter(p0, p1, p2);
        
        let mut corner = Corner::new(t, cc);
        
        // Mark corners outside bounds as border
        if !bounds.contains(&cc) {
            corner.border = true;
        }
        
        mesh.corners.push(corner);
    }
    mesh.num_solid_corners = mesh.corners.len();

    // Step 3: Build point-to-halfedge index for efficient traversal
    let point_to_halfedge = build_point_to_halfedge_index(&triangulation, points.len());

    // Step 4: Build edges and connectivity
    let edge_map = build_edges(&mut mesh, &triangulation, points, bounds);

    // Step 5: Build center-corner connectivity (which corners touch which centers)
    build_center_corner_connectivity(&mut mesh, &triangulation, &point_to_halfedge);

    // Step 6: Build corner-corner adjacency (Voronoi edges connect adjacent corners)
    build_corner_adjacency(&mut mesh, &edge_map);

    // Step 7: Add ghost elements for boundary handling
    add_ghost_elements(&mut mesh, &triangulation, points, bounds);

    // Validate the mesh in debug builds
    #[cfg(debug_assertions)]
    if let Err(e) = mesh.validate() {
        eprintln!("Mesh validation error: {}", e);
    }

    mesh
}

/// Build a mapping from point index to one of its incoming half-edges.
fn build_point_to_halfedge_index(
    triangulation: &delaunator::Triangulation,
    num_points: usize,
) -> Vec<usize> {
    let mut index = vec![delaunator::EMPTY; num_points];
    
    for e in 0..triangulation.triangles.len() {
        let endpoint = triangulation.triangles[next_halfedge(e)];
        // Prefer boundary edges (where halfedges[e] == EMPTY) for boundary points
        if index[endpoint] == delaunator::EMPTY || triangulation.halfedges[e] == delaunator::EMPTY {
            index[endpoint] = e;
        }
    }
    
    index
}

/// Build edges from the triangulation.
/// Returns a map from (min_point, max_point) to edge index.
fn build_edges(
    mesh: &mut DualMesh,
    triangulation: &delaunator::Triangulation,
    points: &[Point],
    _bounds: &BoundingBox,
) -> HashMap<(usize, usize), usize> {
    let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
    let mut processed_halfedges: HashSet<usize> = HashSet::new();

    for e in 0..triangulation.triangles.len() {
        if processed_halfedges.contains(&e) {
            continue;
        }

        let opposite = triangulation.halfedges[e];
        
        // Mark both halfedges as processed
        processed_halfedges.insert(e);
        if opposite != delaunator::EMPTY {
            processed_halfedges.insert(opposite);
        }

        // Get the two endpoint indices (Delaunay edge connects two points)
        let p0 = triangulation.triangles[e];
        let p1 = triangulation.triangles[next_halfedge(e)];

        // Create edge
        let edge_index = mesh.edges.len();
        let mut edge = Edge::new(edge_index);

        // d0, d1 = the two Centers (original points) on either side
        edge.d0 = p0;
        edge.d1 = p1;

        // v0, v1 = the two Corners (triangle circumcenters) at the endpoints
        let t0 = e / 3;  // Triangle containing this halfedge
        edge.v0 = t0;

        if opposite != delaunator::EMPTY {
            let t1 = opposite / 3;  // Triangle on the other side
            edge.v1 = t1;
        } else {
            // Boundary edge - no opposite triangle
            edge.v1 = NONE;
        }

        // Compute midpoint
        let pos0 = points[p0];
        let pos1 = points[p1];
        edge.midpoint = pos0.lerp(&pos1, 0.5);

        // Add edge to centers' border lists
        mesh.centers[p0].borders.push(edge_index);
        mesh.centers[p1].borders.push(edge_index);

        // Add neighbor relationship between centers
        if !mesh.centers[p0].neighbors.contains(&p1) {
            mesh.centers[p0].neighbors.push(p1);
        }
        if !mesh.centers[p1].neighbors.contains(&p0) {
            mesh.centers[p1].neighbors.push(p0);
        }

        // Add edge to corners' protrudes lists
        if edge.v0 != NONE && edge.v0 < mesh.corners.len() {
            mesh.corners[edge.v0].protrudes.push(edge_index);
        }
        if edge.v1 != NONE && edge.v1 < mesh.corners.len() {
            mesh.corners[edge.v1].protrudes.push(edge_index);
        }

        mesh.edges.push(edge);

        // Store in map for later lookup
        let key = if p0 < p1 { (p0, p1) } else { (p1, p0) };
        edge_map.insert(key, edge_index);
    }

    mesh.num_solid_edges = mesh.edges.len();
    edge_map
}

/// Build the connectivity between centers and corners.
/// Each center knows which corners are vertices of its Voronoi polygon.
/// Each corner knows which centers it touches (is a vertex of).
fn build_center_corner_connectivity(
    mesh: &mut DualMesh,
    triangulation: &delaunator::Triangulation,
    point_to_halfedge: &[usize],
) {
    // For each center (point), find all triangles around it
    for center_idx in 0..mesh.centers.len() {
        let start = point_to_halfedge[center_idx];
        if start == delaunator::EMPTY {
            continue;
        }

        let mut corners = Vec::new();
        let mut incoming = start;

        loop {
            let triangle = incoming / 3;
            if triangle < mesh.corners.len() && !corners.contains(&triangle) {
                corners.push(triangle);
                
                // Also update corner's touches list
                if !mesh.corners[triangle].touches.contains(&center_idx) {
                    mesh.corners[triangle].touches.push(center_idx);
                }
            }

            let outgoing = next_halfedge(incoming);
            incoming = triangulation.halfedges[outgoing];

            if incoming == delaunator::EMPTY || incoming == start {
                break;
            }
        }

        mesh.centers[center_idx].corners = corners;
    }
}

/// Build corner-corner adjacency (which corners are connected by Voronoi edges).
fn build_corner_adjacency(
    mesh: &mut DualMesh,
    _edge_map: &HashMap<(usize, usize), usize>,
) {
    // Two corners are adjacent if they share an edge
    for edge in &mesh.edges {
        let v0 = edge.v0;
        let v1 = edge.v1;

        if v0 != NONE && v1 != NONE && v0 < mesh.corners.len() && v1 < mesh.corners.len() {
            if !mesh.corners[v0].adjacent.contains(&v1) {
                mesh.corners[v0].adjacent.push(v1);
            }
            if !mesh.corners[v1].adjacent.contains(&v0) {
                mesh.corners[v1].adjacent.push(v0);
            }
        }
    }
}

/// Add ghost elements to handle boundary cases.
/// Ghost elements ensure that every edge has valid d0/d1 and v0/v1.
fn add_ghost_elements(
    mesh: &mut DualMesh,
    triangulation: &delaunator::Triangulation,
    points: &[Point],
    bounds: &BoundingBox,
) {
    // Create a single ghost center representing "outside the map"
    let ghost_center_idx = mesh.centers.len();
    let ghost_center = Center {
        index: ghost_center_idx,
        position: Point::new(bounds.center().x, bounds.center().y),
        water: true,
        ocean: true,
        ..Default::default()
    };
    mesh.centers.push(ghost_center);

    // Find all boundary edges (where halfedges[e] == EMPTY)
    for e in 0..triangulation.triangles.len() {
        if triangulation.halfedges[e] != delaunator::EMPTY {
            continue;
        }

        // This is a boundary halfedge
        let p0 = triangulation.triangles[e];
        let p1 = triangulation.triangles[next_halfedge(e)];

        // Create a ghost corner for the boundary
        let ghost_corner_idx = mesh.corners.len();
        
        // Position the ghost corner on the boundary edge
        // Project the edge midpoint to the nearest boundary
        let pos0 = points[p0];
        let pos1 = points[p1];
        let midpoint = pos0.lerp(&pos1, 0.5);
        
        // Find intersection with boundary - clamp to bounds
        let ghost_pos = clamp_to_boundary(&midpoint, bounds);

        let mut ghost_corner = Corner::new(ghost_corner_idx, ghost_pos);
        ghost_corner.border = true;
        ghost_corner.water = true;
        ghost_corner.ocean = true;
        mesh.corners.push(ghost_corner);

        // Update the existing edge to reference the ghost corner
        for edge in mesh.edges.iter_mut() {
            if edge.v1 == NONE {
                let (d0, d1) = (edge.d0, edge.d1);
                if (d0 == p0 && d1 == p1) || (d0 == p1 && d1 == p0) {
                    edge.v1 = ghost_corner_idx;
                    
                    // Update adjacency
                    if edge.v0 < mesh.num_solid_corners {
                        mesh.corners[edge.v0].adjacent.push(ghost_corner_idx);
                    }
                    break;
                }
            }
        }

        // Link ghost center with boundary centers
        if !mesh.centers[p0].neighbors.contains(&ghost_center_idx) {
            mesh.centers[p0].neighbors.push(ghost_center_idx);
        }
        if !mesh.centers[p1].neighbors.contains(&ghost_center_idx) {
            mesh.centers[p1].neighbors.push(ghost_center_idx);
        }
        mesh.centers[ghost_center_idx].neighbors.push(p0);
        mesh.centers[ghost_center_idx].neighbors.push(p1);
    }

    // Mark boundary centers
    for &hull_point in &triangulation.hull {
        if hull_point < mesh.num_solid_centers {
            // This center is on the convex hull (boundary)
            mesh.centers[hull_point].ocean = true;
            mesh.centers[hull_point].water = true;
        }
    }
}

// === Halfedge navigation helpers ===

#[inline]
fn next_halfedge(e: usize) -> usize {
    if e % 3 == 2 { e - 2 } else { e + 1 }
}

#[inline]
#[allow(dead_code)]
fn prev_halfedge(e: usize) -> usize {
    if e % 3 == 0 { e + 2 } else { e - 1 }
}

/// Clamp a point to the boundary box, projecting it to the nearest edge.
fn clamp_to_boundary(point: &Point, bounds: &BoundingBox) -> Point {
    let center = bounds.center();
    let dir = *point - center;
    
    if dir.x.abs() < 1e-10 && dir.y.abs() < 1e-10 {
        // Point is at center, push to right edge
        return Point::new(bounds.max.x, center.y);
    }

    // Find where ray from center through point intersects boundary
    let mut t_min = f64::INFINITY;
    
    // Check intersection with each boundary edge
    if dir.x > 0.0 {
        let t = (bounds.max.x - center.x) / dir.x;
        t_min = t_min.min(t);
    } else if dir.x < 0.0 {
        let t = (bounds.min.x - center.x) / dir.x;
        t_min = t_min.min(t);
    }
    
    if dir.y > 0.0 {
        let t = (bounds.max.y - center.y) / dir.y;
        t_min = t_min.min(t);
    } else if dir.y < 0.0 {
        let t = (bounds.min.y - center.y) / dir.y;
        t_min = t_min.min(t);
    }

    // Return intersection point
    Point::new(
        center.x + dir.x * t_min,
        center.y + dir.y * t_min,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pointgen::{PointGeneratorConfig, generate_relaxed_points};

    #[test]
    fn test_build_dual_mesh() {
        let config = PointGeneratorConfig {
            num_points: 50,
            lloyd_iterations: 2,
            ..Default::default()
        };
        let points = generate_relaxed_points(&config);
        let mesh = build_dual_mesh(&points, &config.bounds);

        // Should have same number of centers as input points (plus ghost)
        assert!(mesh.centers.len() >= points.len());
        
        // Should have at least some corners
        assert!(!mesh.corners.is_empty());
        
        // Should have edges
        assert!(!mesh.edges.is_empty());

        // Validate mesh
        assert!(mesh.validate().is_ok());
    }

    #[test]
    fn test_center_corner_connectivity() {
        let config = PointGeneratorConfig {
            num_points: 20,
            lloyd_iterations: 1,
            ..Default::default()
        };
        let points = generate_relaxed_points(&config);
        let mesh = build_dual_mesh(&points, &config.bounds);

        // Each interior center should have at least 3 corners
        for center in mesh.solid_centers() {
            if !center.corners.is_empty() {
                // Non-boundary centers should form a polygon
                assert!(center.corners.len() >= 3 || mesh.is_ghost_center(center.index));
            }
        }
    }
}
