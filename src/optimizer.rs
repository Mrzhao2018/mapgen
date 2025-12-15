//! Automatic parameter optimization for map generation.
//!
//! Uses a hill-climbing algorithm with simulated annealing to find optimal parameters
//! for terrain distribution, biome balance, and coastline quality.
//!
//! The optimizer works directly with GenerationConfig, allowing real parameter testing.

use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use crate::{GenerationConfig, generate_map_with_config, generate_map_stats, MapStats};

/// Target terrain distribution goals.
#[derive(Debug, Clone)]
pub struct TerrainGoals {
    /// Target plains percentage (elevation < 0.3)
    pub plains_target: f64,
    /// Target hills percentage (0.3-0.6)
    pub hills_target: f64,
    /// Target mountains percentage (> 0.6)
    pub mountains_target: f64,
    /// Target land percentage
    pub land_target: f64,
    /// Tolerance for each target (how close is "good enough")
    pub tolerance: f64,
}

impl Default for TerrainGoals {
    fn default() -> Self {
        Self {
            plains_target: 45.0,
            hills_target: 35.0,
            mountains_target: 20.0,
            land_target: 55.0,
            tolerance: 5.0,
        }
    }
}

/// Parameters being optimized.
#[derive(Debug, Clone)]
pub struct OptimizableParams {
    pub redistribution_exponent: f64,
    pub bfs_weight: f64,
    pub noise_frequency: f64,
    pub island_factor: f64,
    pub warp_strength: f64,
    pub shape_octaves: u32,
}

impl Default for OptimizableParams {
    fn default() -> Self {
        Self {
            redistribution_exponent: 0.835,
            bfs_weight: 0.45,
            noise_frequency: 0.0028,
            island_factor: 0.85,  // Larger for more land
            warp_strength: 225.0,
            shape_octaves: 5,
        }
    }
}

impl OptimizableParams {
    /// Create from a GenerationConfig.
    pub fn from_config(config: &GenerationConfig) -> Self {
        Self {
            redistribution_exponent: config.redistribution_exponent,
            bfs_weight: config.bfs_weight,
            noise_frequency: config.elevation_noise_frequency,
            island_factor: config.island_factor,
            warp_strength: config.warp_strength,
            shape_octaves: config.shape_octaves,
        }
    }
    
    /// Apply these parameters to a GenerationConfig.
    pub fn apply_to_config(&self, config: &mut GenerationConfig) {
        config.redistribution_exponent = self.redistribution_exponent;
        config.bfs_weight = self.bfs_weight;
        config.elevation_noise_frequency = self.noise_frequency;
        config.island_factor = self.island_factor;
        config.warp_strength = self.warp_strength;
        config.shape_octaves = self.shape_octaves;
    }
    
    /// Create a GenerationConfig for testing with these parameters.
    pub fn to_test_config(&self, seed: u64) -> GenerationConfig {
        let mut config = GenerationConfig::for_testing(seed);
        self.apply_to_config(&mut config);
        config
    }

    /// Apply small random perturbation to parameters.
    pub fn perturb(&self, rng: &mut impl Rng, step_size: f64) -> Self {
        Self {
            redistribution_exponent: (self.redistribution_exponent + rng.gen_range(-0.15..0.15) * step_size)
                .clamp(0.4, 2.0),
            bfs_weight: (self.bfs_weight + rng.gen_range(-0.1..0.1) * step_size)
                .clamp(0.3, 0.9),
            noise_frequency: (self.noise_frequency + rng.gen_range(-0.0005..0.0005) * step_size)
                .clamp(0.0005, 0.005),
            island_factor: (self.island_factor + rng.gen_range(-0.08..0.08) * step_size)
                .clamp(0.6, 0.95),  // Adjusted range for better land coverage
            warp_strength: (self.warp_strength + rng.gen_range(-40.0..40.0) * step_size)
                .clamp(80.0, 300.0),
            shape_octaves: ((self.shape_octaves as i32 + rng.gen_range(-1..2)) as u32)
                .clamp(3, 7),
        }
    }
}

/// Calculate fitness score (lower is better).
pub fn evaluate_fitness(stats: &MapStats, goals: &TerrainGoals) -> f64 {
    let plains_error = (stats.elevation_stats.plains_percentage - goals.plains_target).abs();
    let hills_error = (stats.elevation_stats.hills_percentage - goals.hills_target).abs();
    let mountains_error = (stats.elevation_stats.mountains_percentage - goals.mountains_target).abs();
    let land_error = (stats.land_percentage - goals.land_target).abs();
    
    // Weighted sum of errors - land coverage is important!
    let terrain_score = plains_error * 1.0 + hills_error * 1.2 + mountains_error * 1.5;
    let land_score = land_error * 1.5;  // Increased weight for land coverage
    
    // Penalize extreme biome distributions
    let biome_penalty = calculate_biome_penalty(stats);
    
    terrain_score + land_score + biome_penalty
}

/// Calculate penalty for poor biome distribution.
fn calculate_biome_penalty(stats: &MapStats) -> f64 {
    use crate::mesh::Biome;
    
    let total = stats.total_regions as f64;
    let mut penalty = 0.0;
    
    // Penalize if any single biome (except Ocean) dominates too much
    for (biome, &count) in &stats.biome_counts {
        if *biome == Biome::Ocean {
            continue;
        }
        let pct = (count as f64 / total) * 100.0;
        if pct > 25.0 {
            penalty += (pct - 25.0) * 0.5;
        }
    }
    
    // Penalize if too few diverse biomes
    let land_biomes = stats.biome_counts.iter()
        .filter(|(b, &c)| **b != Biome::Ocean && **b != Biome::Lake && c > 10)
        .count();
    if land_biomes < 6 {
        penalty += (6 - land_biomes) as f64 * 5.0;
    }
    
    penalty
}

/// Result of optimization.
#[derive(Debug, Clone)]
pub struct OptimizationResult {
    pub best_params: OptimizableParams,
    pub best_fitness: f64,
    pub iterations: usize,
    pub final_stats: Option<MapStats>,
}

/// Run optimization to find best parameters.
pub fn optimize_parameters(
    goals: &TerrainGoals,
    max_iterations: usize,
    samples_per_iteration: usize,
) -> OptimizationResult {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut best_params = OptimizableParams::default();
    let mut best_fitness = f64::INFINITY;
    let mut best_stats = None;
    
    println!("\n=== Starting Parameter Optimization ===");
    println!("Goals: Plains {:.0}%, Hills {:.0}%, Mountains {:.0}%, Land {:.0}%",
        goals.plains_target, goals.hills_target, goals.mountains_target, goals.land_target);
    println!("Max iterations: {}, Samples per iteration: {}\n", max_iterations, samples_per_iteration);
    
    // Start with multiple random starting points
    let mut current_params = best_params.clone();
    let mut step_size = 1.0;
    
    for iter in 0..max_iterations {
        // Evaluate current parameters with multiple samples
        let mut total_fitness = 0.0;
        let mut sample_stats = None;
        
        for sample in 0..samples_per_iteration {
            let seed = 10000 + iter as u64 * 100 + sample as u64;
            let stats = evaluate_params(&current_params, seed);
            let fitness = evaluate_fitness(&stats, goals);
            total_fitness += fitness;
            
            if sample == 0 {
                sample_stats = Some(stats);
            }
        }
        
        let avg_fitness = total_fitness / samples_per_iteration as f64;
        
        // Update best if improved
        if avg_fitness < best_fitness {
            best_fitness = avg_fitness;
            best_params = current_params.clone();
            best_stats = sample_stats;
            
            println!("Iter {}: New best fitness = {:.2}", iter, best_fitness);
            println!("  Params: exp={:.2}, bfs={:.2}, noise={:.4}, island={:.2}, warp={:.0}, oct={}",
                best_params.redistribution_exponent,
                best_params.bfs_weight,
                best_params.noise_frequency,
                best_params.island_factor,
                best_params.warp_strength,
                best_params.shape_octaves);
            
            if let Some(ref s) = best_stats {
                println!("  Terrain: Plains {:.1}%, Hills {:.1}%, Mountains {:.1}%",
                    s.elevation_stats.plains_percentage,
                    s.elevation_stats.hills_percentage,
                    s.elevation_stats.mountains_percentage);
            }
        }
        
        // Generate new candidate by perturbing best
        current_params = best_params.perturb(&mut rng, step_size);
        
        // Decrease step size over time (simulated annealing)
        step_size *= 0.98;
        step_size = step_size.max(0.3);
        
        // Early termination if good enough
        if best_fitness < goals.tolerance * 2.5 {
            println!("\nEarly termination: fitness {:.2} within tolerance", best_fitness);
            break;
        }
    }
    
    println!("\n=== Optimization Complete ===");
    println!("Best fitness: {:.2}", best_fitness);
    println!("Best parameters:");
    println!("  redistribution_exponent: {:.3}", best_params.redistribution_exponent);
    println!("  bfs_weight: {:.3}", best_params.bfs_weight);
    println!("  noise_frequency: {:.5}", best_params.noise_frequency);
    println!("  island_factor: {:.3}", best_params.island_factor);
    println!("  warp_strength: {:.1}", best_params.warp_strength);
    println!("  shape_octaves: {}", best_params.shape_octaves);
    
    OptimizationResult {
        best_params,
        best_fitness,
        iterations: max_iterations,
        final_stats: best_stats,
    }
}

/// Evaluate parameters by generating a map and collecting stats.
/// Now properly uses GenerationConfig to pass all parameters.
fn evaluate_params(params: &OptimizableParams, seed: u64) -> MapStats {
    // Create test config with the parameters applied
    let config = params.to_test_config(seed);
    
    // Generate map with full parameter propagation
    let mesh = generate_map_with_config(&config);
    generate_map_stats(&mesh)
}

/// Quick search to find a good starting point.
pub fn quick_search(goals: &TerrainGoals) -> OptimizableParams {
    println!("\n=== Quick Parameter Search ===");
    
    // Grid search over key parameters
    let exponents = [0.7, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5];
    let bfs_weights = [0.6, 0.7, 0.8];
    
    let mut best_params = OptimizableParams::default();
    let mut best_fitness = f64::INFINITY;
    
    for &exp in &exponents {
        for &bfs in &bfs_weights {
            let params = OptimizableParams {
                redistribution_exponent: exp,
                bfs_weight: bfs,
                ..Default::default()
            };
            
            let stats = evaluate_params(&params, 12345);
            let fitness = evaluate_fitness(&stats, goals);
            
            println!("exp={:.2}, bfs={:.2} -> fitness={:.2} (P:{:.1}%, H:{:.1}%, M:{:.1}%)",
                exp, bfs, fitness,
                stats.elevation_stats.plains_percentage,
                stats.elevation_stats.hills_percentage,
                stats.elevation_stats.mountains_percentage);
            
            if fitness < best_fitness {
                best_fitness = fitness;
                best_params = params;
            }
        }
    }
    
    println!("\nBest quick search: exp={:.2}, bfs={:.2}, fitness={:.2}",
        best_params.redistribution_exponent, best_params.bfs_weight, best_fitness);
    
    best_params
}

/// Print recommended changes - now outputs GenerationConfig format.
pub fn print_recommendations(result: &OptimizationResult) {
    println!("\n╔══════════════════════════════════════════════════════════════╗");
    println!("║                  Optimization Results                        ║");
    println!("╚══════════════════════════════════════════════════════════════╝\n");
    
    println!("Best Parameters (GenerationConfig format):");
    println!("─────────────────────────────────────────────");
    println!("  redistribution_exponent: {:.3}", result.best_params.redistribution_exponent);
    println!("  bfs_weight:              {:.3}", result.best_params.bfs_weight);
    println!("  elevation_noise_freq:    {:.5}", result.best_params.noise_frequency);
    println!("  island_factor:           {:.3}", result.best_params.island_factor);
    println!("  warp_strength:           {:.1}", result.best_params.warp_strength);
    println!("  shape_octaves:           {}", result.best_params.shape_octaves);
    
    if let Some(ref stats) = result.final_stats {
        println!("\nResulting Terrain Distribution:");
        println!("─────────────────────────────────────────────");
        println!("  Plains (<0.3):     {:.1}%", stats.elevation_stats.plains_percentage);
        println!("  Hills (0.3-0.6):   {:.1}%", stats.elevation_stats.hills_percentage);
        println!("  Mountains (>0.6):  {:.1}%", stats.elevation_stats.mountains_percentage);
        println!("  Land Coverage:     {:.1}%", stats.land_percentage);
    }
    
    println!("\nTo use these parameters in code:");
    println!("─────────────────────────────────────────────");
    println!("```rust");
    println!("let config = GenerationConfig {{");
    println!("    redistribution_exponent: {:.3},", result.best_params.redistribution_exponent);
    println!("    bfs_weight: {:.3},", result.best_params.bfs_weight);
    println!("    elevation_noise_frequency: {:.5},", result.best_params.noise_frequency);
    println!("    island_factor: {:.3},", result.best_params.island_factor);
    println!("    warp_strength: {:.1},", result.best_params.warp_strength);
    println!("    shape_octaves: {},", result.best_params.shape_octaves);
    println!("    ..Default::default()");
    println!("}};");
    println!("let mesh = generate_map_with_config(&config);");
    println!("```");
}

/// Run a full optimization and return the best config.
pub fn optimize_and_create_config(goals: &TerrainGoals, seed: u64) -> GenerationConfig {
    let result = optimize_parameters(goals, 60, 3);
    let mut config = GenerationConfig::with_seed(seed);
    result.best_params.apply_to_config(&mut config);
    config
}
