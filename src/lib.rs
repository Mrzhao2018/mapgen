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
pub mod optimizer;
pub mod config;

use geometry::BoundingBox;
use pointgen::PointGeneratorConfig;
use graph::build_dual_mesh;
use island::{IslandConfig, IslandShape, generate_island};
use elevation::{ElevationConfig, generate_elevation};
use rivers::{RiverConfig, generate_rivers, calculate_moisture};
use biomes::assign_biomes;
use mesh::DualMesh;

// Re-export the unified config
pub use config::GenerationConfig;

/// Legacy map configuration (deprecated, use GenerationConfig instead).
#[derive(Debug, Clone)]
#[deprecated(note = "Use GenerationConfig instead")]
pub struct MapConfig {
    pub bounds: BoundingBox,
    pub num_points: usize,
    pub seed: u64,
    pub lloyd_iterations: u32,
    pub island_shape: IslandShape,
    pub island_factor: f64,
    pub num_rivers: usize,
}

#[allow(deprecated)]
impl Default for MapConfig {
    fn default() -> Self {
        Self {
            bounds: BoundingBox::new(0.0, 0.0, 2000.0, 2000.0),
            num_points: 4000,
            seed: 12345,
            lloyd_iterations: 2,
            island_shape: IslandShape::DomainWarp,
            island_factor: 0.65,
            num_rivers: 80,
        }
    }
}

/// Generate a complete map using the unified configuration.
pub fn generate_map_with_config(config: &GenerationConfig) -> DualMesh {
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

    // Step 3: Generate island shape with full config
    let island_config = IslandConfig {
        seed: config.seed as u32,
        noise_frequency: config.coastline_noise_frequency,
        noise_amplitude: config.coastline_noise_amplitude,
        island_factor: config.island_factor,
        shape: config.island_shape,
    };
    generate_island(&mut mesh, &island_config, &config.bounds);
    
    // Temporarily disabled to debug the "cut" issue
    // Remove small islands
    // island::remove_small_islands(&mut mesh, config.min_island_cells, config.min_island_fraction);

    // Step 4: Calculate elevation with unified config
    let elevation_config = ElevationConfig {
        redistribution_exponent: config.redistribution_exponent,
        lake_elevation_factor: config.lake_elevation_factor,
        use_hybrid: true,
        bfs_weight: config.bfs_weight,
        noise_frequency: config.elevation_noise_frequency,
        fill_pits: config.fill_pits,
        seed: config.seed as u32,
        use_ridged: config.use_ridged_noise,
    };
    generate_elevation(&mut mesh, &elevation_config);

    // Step 5: Generate rivers
    let river_config = RiverConfig {
        seed: config.seed,
        num_rivers: config.num_rivers,
        min_source_elevation: config.river_min_elevation,
        max_source_elevation: config.river_max_elevation,
    };
    generate_rivers(&mut mesh, &river_config);

    // Step 6: Calculate moisture
    calculate_moisture(&mut mesh);

    // Step 7: Assign biomes
    assign_biomes(&mut mesh);

    mesh
}

/// Generate a complete map (legacy API, converts to new config internally).
#[allow(deprecated)]
pub fn generate_map(config: &MapConfig) -> DualMesh {
    // Convert legacy config to new unified config
    let unified = GenerationConfig {
        bounds: config.bounds,
        num_points: config.num_points,
        seed: config.seed,
        lloyd_iterations: config.lloyd_iterations,
        island_shape: config.island_shape,
        island_factor: config.island_factor,
        num_rivers: config.num_rivers,
        ..Default::default()
    };
    
    generate_map_with_config(&unified)
}

/// Statistics report for generated map.
#[derive(Debug, Clone)]
pub struct MapStats {
    pub total_regions: usize,
    pub land_regions: usize,
    pub ocean_regions: usize,
    pub lake_regions: usize,
    pub coast_regions: usize,
    pub land_percentage: f64,
    pub total_corners: usize,
    pub total_edges: usize,
    pub river_edges: usize,
    pub biome_counts: std::collections::HashMap<mesh::Biome, usize>,
    pub elevation_stats: ElevationStats,
}

#[derive(Debug, Clone)]
pub struct ElevationStats {
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    pub plains_percentage: f64,    // elevation < 0.3
    pub hills_percentage: f64,     // 0.3 <= elevation < 0.6
    pub mountains_percentage: f64, // elevation >= 0.6
}

/// Generate statistics report for a mesh.
pub fn generate_map_stats(mesh: &DualMesh) -> MapStats {
    use std::collections::HashMap;
    
    let total_regions = mesh.num_solid_centers;
    let land_regions = mesh.solid_centers().filter(|c| !c.water).count();
    let ocean_regions = mesh.solid_centers().filter(|c| c.ocean).count();
    let lake_regions = mesh.solid_centers().filter(|c| c.water && !c.ocean).count();
    let coast_regions = mesh.solid_centers().filter(|c| c.coast).count();
    let land_percentage = (land_regions as f64 / total_regions as f64) * 100.0;
    
    let river_edges = mesh.edges.iter().filter(|e| e.river > 0).count();
    
    // Biome counts
    let mut biome_counts: HashMap<mesh::Biome, usize> = HashMap::new();
    for center in mesh.solid_centers() {
        *biome_counts.entry(center.biome).or_insert(0) += 1;
    }
    
    // Elevation stats for land only
    let land_elevations: Vec<f64> = mesh.solid_centers()
        .filter(|c| !c.water)
        .map(|c| c.elevation)
        .collect();
    
    let elevation_stats = if !land_elevations.is_empty() {
        let min = land_elevations.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = land_elevations.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean = land_elevations.iter().sum::<f64>() / land_elevations.len() as f64;
        
        let plains = land_elevations.iter().filter(|&&e| e < 0.3).count();
        let hills = land_elevations.iter().filter(|&&e| e >= 0.3 && e < 0.6).count();
        let mountains = land_elevations.iter().filter(|&&e| e >= 0.6).count();
        let total = land_elevations.len() as f64;
        
        ElevationStats {
            min,
            max,
            mean,
            plains_percentage: (plains as f64 / total) * 100.0,
            hills_percentage: (hills as f64 / total) * 100.0,
            mountains_percentage: (mountains as f64 / total) * 100.0,
        }
    } else {
        ElevationStats {
            min: 0.0, max: 0.0, mean: 0.0,
            plains_percentage: 0.0, hills_percentage: 0.0, mountains_percentage: 0.0,
        }
    };
    
    MapStats {
        total_regions,
        land_regions,
        ocean_regions,
        lake_regions,
        coast_regions,
        land_percentage,
        total_corners: mesh.num_solid_corners,
        total_edges: mesh.edges.len(),
        river_edges,
        biome_counts,
        elevation_stats,
    }
}

/// Print map statistics report.
pub fn print_map_report(stats: &MapStats) {
    println!("\n=== Map Generation Report ===");
    println!("Regions: {} total", stats.total_regions);
    println!("  - Land: {} ({:.1}%)", stats.land_regions, stats.land_percentage);
    println!("  - Ocean: {}", stats.ocean_regions);
    println!("  - Lakes: {}", stats.lake_regions);
    println!("  - Coast: {}", stats.coast_regions);
    println!("Corners: {}", stats.total_corners);
    println!("Edges: {} (rivers: {})", stats.total_edges, stats.river_edges);
    
    println!("\nElevation Distribution (land only):");
    println!("  - Range: {:.3} to {:.3} (mean: {:.3})", 
        stats.elevation_stats.min, stats.elevation_stats.max, stats.elevation_stats.mean);
    println!("  - Plains (<0.3): {:.1}%", stats.elevation_stats.plains_percentage);
    println!("  - Hills (0.3-0.6): {:.1}%", stats.elevation_stats.hills_percentage);
    println!("  - Mountains (>0.6): {:.1}%", stats.elevation_stats.mountains_percentage);
    
    println!("\nBiome Distribution:");
    let mut biomes: Vec<_> = stats.biome_counts.iter().collect();
    biomes.sort_by(|a, b| b.1.cmp(a.1));
    for (biome, count) in biomes {
        let pct = (*count as f64 / stats.total_regions as f64) * 100.0;
        println!("  - {:?}: {} ({:.1}%)", biome, count, pct);
    }
    println!("=============================\n");
}
