//! Unified configuration system for map generation.
//!
//! All generation parameters are centralized here, allowing easy adjustment
//! and optimization of the entire generation pipeline.

use crate::geometry::BoundingBox;
use crate::island::IslandShape;

/// Complete unified configuration for map generation.
/// Contains all tunable parameters across all generation stages.
#[derive(Debug, Clone)]
pub struct GenerationConfig {
    // ===== Basic Map Settings =====
    /// Bounding box for the map.
    pub bounds: BoundingBox,
    /// Number of Voronoi cells.
    pub num_points: usize,
    /// Random seed for generation.
    pub seed: u64,
    /// Lloyd relaxation iterations.
    pub lloyd_iterations: u32,
    
    // ===== Island Shape Settings =====
    /// Island shape type.
    pub island_shape: IslandShape,
    /// Island size factor (0.0 = tiny, 1.0 = fills map).
    pub island_factor: f64,
    /// Noise frequency for coastline detail.
    pub coastline_noise_frequency: f64,
    /// Noise amplitude for coastline variation.
    pub coastline_noise_amplitude: f64,
    /// Domain warp strength (only for DomainWarp shape).
    pub warp_strength: f64,
    /// Domain warp iterations.
    pub warp_iterations: u32,
    /// Shape noise octaves.
    pub shape_octaves: u32,
    
    // ===== Elevation Settings =====
    /// Redistribution exponent (higher = more lowlands, lower = more mountains).
    /// Uses formula: elevation = x^exponent
    pub redistribution_exponent: f64,
    /// Weight for BFS distance field vs noise (0-1).
    pub bfs_weight: f64,
    /// Noise frequency for elevation modulation.
    pub elevation_noise_frequency: f64,
    /// Use ridged multi-fractal for mountains.
    pub use_ridged_noise: bool,
    /// Enable pit-filling algorithm.
    pub fill_pits: bool,
    /// Lake elevation factor.
    pub lake_elevation_factor: f64,
    
    // ===== River Settings =====
    /// Number of rivers to generate.
    pub num_rivers: usize,
    /// Minimum elevation for river sources.
    pub river_min_elevation: f64,
    /// Maximum elevation for river sources.
    pub river_max_elevation: f64,
    
    // ===== Small Island Removal =====
    /// Minimum cells for an island to survive.
    pub min_island_cells: usize,
    /// Minimum island size as fraction of largest island.
    pub min_island_fraction: f64,
}

impl Default for GenerationConfig {
    fn default() -> Self {
        Self {
            // Basic
            bounds: BoundingBox::new(0.0, 0.0, 2000.0, 2000.0),
            num_points: 4000,
            seed: 12345,
            lloyd_iterations: 2,
            
            // Island shape - larger island for more land coverage
            island_shape: IslandShape::DomainWarp,
            island_factor: 0.85,  // Larger for ~50% land
            coastline_noise_frequency: 0.005,
            coastline_noise_amplitude: 0.5,
            warp_strength: 225.0,
            warp_iterations: 2,
            shape_octaves: 5,
            
            // Elevation - optimized for ~40% plains, ~35% hills, ~25% mountains
            redistribution_exponent: 0.835,
            bfs_weight: 0.45,
            elevation_noise_frequency: 0.0028,
            use_ridged_noise: true,
            fill_pits: true,
            lake_elevation_factor: 0.08,
            
            // Rivers
            num_rivers: 80,
            river_min_elevation: 0.4,
            river_max_elevation: 0.9,
            
            // Island cleanup - more aggressive to remove noise fragments
            min_island_cells: 50,      // Minimum 50 cells to survive
            min_island_fraction: 0.15, // Must be at least 15% of largest island
        }
    }
}

impl GenerationConfig {
    /// Create config with custom seed.
    pub fn with_seed(seed: u64) -> Self {
        Self { seed, ..Default::default() }
    }
    
    /// Create a smaller config for faster testing/optimization.
    pub fn for_testing(seed: u64) -> Self {
        Self {
            bounds: BoundingBox::new(0.0, 0.0, 1500.0, 1500.0),
            num_points: 2000,
            seed,
            num_rivers: 30,
            ..Default::default()
        }
    }
    
    /// Create high-quality config for final rendering.
    pub fn high_quality(seed: u64) -> Self {
        Self {
            bounds: BoundingBox::new(0.0, 0.0, 2000.0, 2000.0),
            num_points: 8000,
            seed,
            lloyd_iterations: 3,
            num_rivers: 100,
            ..Default::default()
        }
    }
    
    /// Apply optimized terrain parameters.
    pub fn with_terrain_params(
        mut self,
        redistribution_exponent: f64,
        bfs_weight: f64,
        elevation_noise_frequency: f64,
    ) -> Self {
        self.redistribution_exponent = redistribution_exponent;
        self.bfs_weight = bfs_weight;
        self.elevation_noise_frequency = elevation_noise_frequency;
        self
    }
    
    /// Apply island shape parameters.
    pub fn with_island_params(
        mut self,
        island_factor: f64,
        warp_strength: f64,
        shape_octaves: u32,
    ) -> Self {
        self.island_factor = island_factor;
        self.warp_strength = warp_strength;
        self.shape_octaves = shape_octaves;
        self
    }
}

/// Presets for common terrain distributions.
pub mod presets {
    use super::GenerationConfig;
    
    /// More plains, gentle terrain.
    pub fn gentle_terrain(seed: u64) -> GenerationConfig {
        GenerationConfig {
            seed,
            redistribution_exponent: 1.2,
            bfs_weight: 0.7,
            ..Default::default()
        }
    }
    
    /// Balanced terrain (~40% plains, ~35% hills, ~25% mountains).
    pub fn balanced_terrain(seed: u64) -> GenerationConfig {
        GenerationConfig {
            seed,
            redistribution_exponent: 0.75,
            bfs_weight: 0.55,
            ..Default::default()
        }
    }
    
    /// More mountains, rugged terrain.
    pub fn mountainous_terrain(seed: u64) -> GenerationConfig {
        GenerationConfig {
            seed,
            redistribution_exponent: 0.5,
            bfs_weight: 0.4,
            ..Default::default()
        }
    }
    
    /// Archipelago style (many small islands).
    pub fn archipelago(seed: u64) -> GenerationConfig {
        GenerationConfig {
            seed,
            island_shape: super::IslandShape::Noise,
            island_factor: 0.5,
            redistribution_exponent: 0.9,
            ..Default::default()
        }
    }
}
