//! Polygonal Map Generator
//!
//! Based on the algorithm by Amit Patel (Red Blob Games):
//! https://www.redblobgames.com/maps/mapgen2/
//!
//! This implementation uses:
//! - Arena-based data structures (no Rc<RefCell<T>>)
//! - Voronoi/Delaunay dual graph system
//! - Gameplay constraints drive generation approach

pub mod geometry;
pub mod mesh;
pub mod pointgen;
pub mod graph;
pub mod island;
pub mod elevation;
pub mod rivers;
pub mod biomes;
pub mod visualizer;

use geometry::BoundingBox;
use pointgen::PointGeneratorConfig;
use graph::build_dual_mesh;
use island::{IslandConfig, IslandShape, generate_island};
use elevation::{ElevationConfig, generate_elevation};
use rivers::{RiverConfig, generate_rivers, calculate_moisture};
use biomes::assign_biomes;
use mesh::DualMesh;

/// Complete map generation configuration.
#[derive(Debug, Clone)]
pub struct MapConfig {
    pub bounds: BoundingBox,
    pub num_points: usize,
    pub seed: u64,
    pub lloyd_iterations: u32,
    pub island_shape: IslandShape,
    pub island_factor: f64,
    pub num_rivers: usize,
}

impl Default for MapConfig {
    fn default() -> Self {
        Self {
            bounds: BoundingBox::new(0.0, 0.0, 1000.0, 1000.0),
            num_points: 1000,
            seed: 12345,
            lloyd_iterations: 2,
            island_shape: IslandShape::Radial,
            island_factor: 0.7,
            num_rivers: 30,
        }
    }
}

/// Generate a complete map.
pub fn generate_map(config: &MapConfig) -> DualMesh {
    // Step 1: Generate and relax points
    let point_config = PointGeneratorConfig {
        bounds: config.bounds,
        num_points: config.num_points,
        seed: config.seed,
        lloyd_iterations: config.lloyd_iterations,
        lloyd_omega: 1.0,
    };
    let points = pointgen::generate_relaxed_points(&point_config);

    // Step 2: Build dual graph
    let mut mesh = build_dual_mesh(&points, &config.bounds);

    // Step 3: Generate island shape
    let island_config = IslandConfig {
        seed: config.seed as u32,
        noise_frequency: 0.005,
        noise_amplitude: 0.5,
        island_factor: config.island_factor,
        shape: config.island_shape,
    };
    generate_island(&mut mesh, &island_config, &config.bounds);

    // Step 4: Calculate elevation
    let elevation_config = ElevationConfig::default();
    generate_elevation(&mut mesh, &elevation_config);

    // Step 5: Generate rivers
    let river_config = RiverConfig {
        seed: config.seed,
        num_rivers: config.num_rivers,
        ..Default::default()
    };
    generate_rivers(&mut mesh, &river_config);

    // Step 6: Calculate moisture
    calculate_moisture(&mut mesh);

    // Step 7: Assign biomes
    assign_biomes(&mut mesh);

    mesh
}
