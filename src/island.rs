//! Island shape generation using noise and radial masks.
//!
//! This module determines which regions are land vs water,
//! and distinguishes ocean (connected to boundary) from lakes.

use std::collections::VecDeque;
use noise::{NoiseFn, Perlin};
use crate::geometry::{Point, BoundingBox};
use crate::mesh::{DualMesh, NONE};

/// Configuration for island shape generation.
#[derive(Debug, Clone)]
pub struct IslandConfig {
    /// Random seed for noise generation.
    pub seed: u32,
    /// Noise frequency (higher = more detailed coastline).
    pub noise_frequency: f64,
    /// Noise amplitude (how much noise affects the shape).
    pub noise_amplitude: f64,
    /// Island size factor (0.0 = tiny island, 1.0 = fills map).
    pub island_factor: f64,
    /// Shape type.
    pub shape: IslandShape,
}

impl Default for IslandConfig {
    fn default() -> Self {
        Self {
            seed: 12345,
            noise_frequency: 0.005,
            noise_amplitude: 0.5,
            island_factor: 0.7,
            shape: IslandShape::Radial,
        }
    }
}

/// Different island shape types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IslandShape {
    /// Circular island with noise.
    Radial,
    /// Square island with noise.
    Square,
    /// Perlin noise only (archipelago-like).
    Noise,
    /// Blob shape (smooth organic).
    Blob,
}

/// Shape function that determines if a point should be land.
/// Returns a value from 0.0 (definitely water) to 1.0 (definitely land).
pub trait ShapeFunction {
    fn evaluate(&self, point: &Point, bounds: &BoundingBox) -> f64;
}

/// Radial island shape with noise perturbation.
pub struct RadialIsland {
    perlin: Perlin,
    config: IslandConfig,
}

impl RadialIsland {
    pub fn new(config: &IslandConfig) -> Self {
        let perlin = Perlin::new(config.seed);
        Self {
            perlin,
            config: config.clone(),
        }
    }
}

impl ShapeFunction for RadialIsland {
    fn evaluate(&self, point: &Point, bounds: &BoundingBox) -> f64 {
        let center = bounds.center();
        let max_radius = bounds.width().min(bounds.height()) / 2.0;

        // Normalize coordinates to [-1, 1]
        let nx = (point.x - center.x) / max_radius;
        let ny = (point.y - center.y) / max_radius;

        // Distance from center (0 at center, 1 at edge)
        let distance = (nx * nx + ny * ny).sqrt();

        // Add noise perturbation
        let noise_val = self.perlin.get([
            point.x * self.config.noise_frequency,
            point.y * self.config.noise_frequency,
        ]);
        let noise_contribution = noise_val * self.config.noise_amplitude;

        // Base shape: radial gradient
        let base = 1.0 - distance / self.config.island_factor;
        
        // Combine base shape with noise
        let value = base + noise_contribution;

        // Clamp to [0, 1]
        value.clamp(0.0, 1.0)
    }
}

/// Square island shape with noise perturbation.
pub struct SquareIsland {
    perlin: Perlin,
    config: IslandConfig,
}

impl SquareIsland {
    pub fn new(config: &IslandConfig) -> Self {
        let perlin = Perlin::new(config.seed);
        Self {
            perlin,
            config: config.clone(),
        }
    }
}

impl ShapeFunction for SquareIsland {
    fn evaluate(&self, point: &Point, bounds: &BoundingBox) -> f64 {
        let center = bounds.center();
        let half_width = bounds.width() / 2.0;
        let half_height = bounds.height() / 2.0;

        // Normalize to [-1, 1]
        let nx = (point.x - center.x) / half_width;
        let ny = (point.y - center.y) / half_height;

        // Use Chebyshev distance (max of abs)
        let distance = nx.abs().max(ny.abs());

        // Add noise
        let noise_val = self.perlin.get([
            point.x * self.config.noise_frequency,
            point.y * self.config.noise_frequency,
        ]);
        let noise_contribution = noise_val * self.config.noise_amplitude;

        let base = 1.0 - distance / self.config.island_factor;
        let value = base + noise_contribution;

        value.clamp(0.0, 1.0)
    }
}

/// Pure noise-based shape (creates archipelago).
pub struct NoiseIsland {
    perlin: Perlin,
    config: IslandConfig,
}

impl NoiseIsland {
    pub fn new(config: &IslandConfig) -> Self {
        let perlin = Perlin::new(config.seed);
        Self {
            perlin,
            config: config.clone(),
        }
    }
}

impl ShapeFunction for NoiseIsland {
    fn evaluate(&self, point: &Point, bounds: &BoundingBox) -> f64 {
        let center = bounds.center();
        let max_dist = bounds.width().min(bounds.height()) / 2.0;

        // Edge falloff
        let nx = (point.x - center.x) / max_dist;
        let ny = (point.y - center.y) / max_dist;
        let edge_distance = (nx * nx + ny * ny).sqrt();
        let edge_factor = 1.0 - edge_distance.powf(2.0);

        // Multi-octave noise
        let noise1 = self.perlin.get([
            point.x * self.config.noise_frequency,
            point.y * self.config.noise_frequency,
        ]);
        let noise2 = self.perlin.get([
            point.x * self.config.noise_frequency * 2.0 + 100.0,
            point.y * self.config.noise_frequency * 2.0 + 100.0,
        ]) * 0.5;
        let noise3 = self.perlin.get([
            point.x * self.config.noise_frequency * 4.0 + 200.0,
            point.y * self.config.noise_frequency * 4.0 + 200.0,
        ]) * 0.25;

        let combined_noise = (noise1 + noise2 + noise3) / 1.75;
        let value = (combined_noise + 0.5) * edge_factor * self.config.island_factor;

        value.clamp(0.0, 1.0)
    }
}

/// Blob-shaped island (smooth, organic).
pub struct BlobIsland {
    perlin: Perlin,
    config: IslandConfig,
    bumps: usize,
    angle_offset: f64,
}

impl BlobIsland {
    pub fn new(config: &IslandConfig) -> Self {
        let perlin = Perlin::new(config.seed);
        // Use seed to determine blob characteristics
        let bumps = 3 + (config.seed % 5) as usize;
        let angle_offset = (config.seed as f64 * 0.1) % (2.0 * std::f64::consts::PI);
        
        Self {
            perlin,
            config: config.clone(),
            bumps,
            angle_offset,
        }
    }
}

impl ShapeFunction for BlobIsland {
    fn evaluate(&self, point: &Point, bounds: &BoundingBox) -> f64 {
        let center = bounds.center();
        let max_radius = bounds.width().min(bounds.height()) / 2.0;

        let dx = point.x - center.x;
        let dy = point.y - center.y;
        let distance = (dx * dx + dy * dy).sqrt();
        let angle = dy.atan2(dx);

        // Blob radius varies with angle
        let blob_radius = max_radius * self.config.island_factor * (
            0.7 + 0.3 * ((self.bumps as f64) * angle + self.angle_offset).sin()
        );

        // Add noise
        let noise_val = self.perlin.get([
            point.x * self.config.noise_frequency,
            point.y * self.config.noise_frequency,
        ]);
        let noise_radius = blob_radius * (1.0 + noise_val * self.config.noise_amplitude * 0.5);

        if distance < noise_radius {
            1.0 - (distance / noise_radius).powf(2.0)
        } else {
            0.0
        }
    }
}

/// Create a shape function based on configuration.
pub fn create_shape_function(config: &IslandConfig) -> Box<dyn ShapeFunction> {
    match config.shape {
        IslandShape::Radial => Box::new(RadialIsland::new(config)),
        IslandShape::Square => Box::new(SquareIsland::new(config)),
        IslandShape::Noise => Box::new(NoiseIsland::new(config)),
        IslandShape::Blob => Box::new(BlobIsland::new(config)),
    }
}

/// Assign water/land to all corners based on the shape function.
pub fn assign_corner_water(
    mesh: &mut DualMesh,
    shape: &dyn ShapeFunction,
    bounds: &BoundingBox,
    water_threshold: f64,
) {
    for corner in mesh.corners.iter_mut() {
        let value = shape.evaluate(&corner.position, bounds);
        corner.water = value < water_threshold;
        
        // Border corners are always water/ocean
        if corner.border {
            corner.water = true;
            corner.ocean = true;
        }
    }
}

/// Assign water/land to centers based on their corners.
/// A center is water if all its corners are water.
/// A center is coast if it has both land and water corners.
pub fn assign_center_water(mesh: &mut DualMesh) {
    // First pass: determine water status
    for i in 0..mesh.centers.len() {
        let corners = mesh.centers[i].corners.clone();
        
        if corners.is_empty() {
            // Ghost or degenerate center
            continue;
        }

        let num_water = corners.iter()
            .filter(|&&c| c < mesh.corners.len() && mesh.corners[c].water)
            .count();
        
        // A center is water if majority of corners are water
        mesh.centers[i].water = num_water > corners.len() / 2;
    }

    // Mark boundary centers as water
    for i in 0..mesh.num_solid_centers {
        let is_boundary = mesh.centers[i].neighbors.iter()
            .any(|&n| n >= mesh.num_solid_centers);
        
        if is_boundary {
            mesh.centers[i].water = true;
            mesh.centers[i].ocean = true;
        }
    }
}

/// Use flood fill to identify ocean (water connected to boundary) vs lakes.
pub fn assign_ocean(mesh: &mut DualMesh) {
    // Start from all boundary corners and flood fill
    let mut queue = VecDeque::new();
    
    // Initialize: all border corners are ocean
    for i in 0..mesh.corners.len() {
        if mesh.corners[i].border && mesh.corners[i].water {
            mesh.corners[i].ocean = true;
            queue.push_back(i);
        }
    }

    // Also start from corners that touch the ghost center
    let ghost_center_idx = mesh.num_solid_centers;
    if ghost_center_idx < mesh.centers.len() {
        for &corner_idx in &mesh.centers[ghost_center_idx].corners.clone() {
            if corner_idx < mesh.corners.len() && mesh.corners[corner_idx].water {
                mesh.corners[corner_idx].ocean = true;
                if !queue.contains(&corner_idx) {
                    queue.push_back(corner_idx);
                }
            }
        }
    }

    // Flood fill through water corners
    while let Some(corner_idx) = queue.pop_front() {
        for &adj in &mesh.corners[corner_idx].adjacent.clone() {
            if adj == NONE || adj >= mesh.corners.len() {
                continue;
            }
            if mesh.corners[adj].water && !mesh.corners[adj].ocean {
                mesh.corners[adj].ocean = true;
                queue.push_back(adj);
            }
        }
    }

    // Propagate ocean status to centers
    for i in 0..mesh.centers.len() {
        if mesh.centers[i].water {
            let corners = &mesh.centers[i].corners;
            let has_ocean_corner = corners.iter()
                .any(|&c| c < mesh.corners.len() && mesh.corners[c].ocean);
            mesh.centers[i].ocean = has_ocean_corner;
        }
    }

    // Identify coast centers (land adjacent to ocean)
    for i in 0..mesh.num_solid_centers {
        if !mesh.centers[i].water {
            let neighbors = mesh.centers[i].neighbors.clone();
            let adjacent_to_ocean = neighbors.iter()
                .any(|&n| n < mesh.centers.len() && mesh.centers[n].ocean);
            mesh.centers[i].coast = adjacent_to_ocean;
        }
    }

    // Identify coast corners (touch both land and ocean centers)
    for i in 0..mesh.corners.len() {
        let touches = &mesh.corners[i].touches;
        let touch_land = touches.iter()
            .any(|&c| c < mesh.centers.len() && !mesh.centers[c].water);
        let touch_ocean = touches.iter()
            .any(|&c| c < mesh.centers.len() && mesh.centers[c].ocean);
        mesh.corners[i].coast = touch_land && touch_ocean;
    }
}

/// Generate island shape for the mesh.
pub fn generate_island(
    mesh: &mut DualMesh,
    config: &IslandConfig,
    bounds: &BoundingBox,
) {
    let shape = create_shape_function(config);
    
    // Water threshold: values below this are water
    let water_threshold = 0.3;
    
    assign_corner_water(mesh, shape.as_ref(), bounds, water_threshold);
    assign_center_water(mesh);
    assign_ocean(mesh);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pointgen::{PointGeneratorConfig, generate_relaxed_points};
    use crate::graph::build_dual_mesh;

    #[test]
    fn test_island_generation() {
        let point_config = PointGeneratorConfig {
            num_points: 100,
            lloyd_iterations: 2,
            ..Default::default()
        };
        let points = generate_relaxed_points(&point_config);
        let mut mesh = build_dual_mesh(&points, &point_config.bounds);

        let island_config = IslandConfig {
            shape: IslandShape::Radial,
            island_factor: 0.7,
            ..Default::default()
        };
        generate_island(&mut mesh, &island_config, &point_config.bounds);

        // Should have some land and some water
        let land_count = mesh.solid_centers()
            .filter(|c| !c.water)
            .count();
        let water_count = mesh.solid_centers()
            .filter(|c| c.water)
            .count();

        assert!(land_count > 0, "Should have some land");
        assert!(water_count > 0, "Should have some water");
    }

    #[test]
    fn test_ocean_flood_fill() {
        let point_config = PointGeneratorConfig {
            num_points: 50,
            lloyd_iterations: 1,
            ..Default::default()
        };
        let points = generate_relaxed_points(&point_config);
        let mut mesh = build_dual_mesh(&points, &point_config.bounds);

        let island_config = IslandConfig::default();
        generate_island(&mut mesh, &island_config, &point_config.bounds);

        // All boundary centers should be ocean
        for i in 0..mesh.num_solid_centers {
            if mesh.centers[i].neighbors.iter().any(|&n| n >= mesh.num_solid_centers) {
                assert!(mesh.centers[i].ocean, "Boundary center {} should be ocean", i);
            }
        }
    }
}
