//! Parameter optimization tool for map generation.
//!
//! Run with: cargo run --bin optimize
//!
//! This tool searches for optimal terrain generation parameters
//! to achieve a target terrain distribution.

use mapgen::optimizer::{TerrainGoals, quick_search, optimize_parameters, print_recommendations};
use mapgen::{GenerationConfig, generate_map_with_config, generate_map_stats, print_map_report};

fn main() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║        Polygonal Map Generator - Parameter Optimizer         ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    
    // Define target terrain distribution
    let goals = TerrainGoals {
        plains_target: 40.0,      // 40% plains (elevation < 0.3)
        hills_target: 35.0,       // 35% hills (0.3-0.6)
        mountains_target: 25.0,   // 25% mountains (> 0.6)
        land_target: 50.0,        // 50% land coverage
        tolerance: 5.0,           // ±5% tolerance
    };
    
    println!("\n┌─────────────────────────────────────────────┐");
    println!("│ Target Distribution                         │");
    println!("├─────────────────────────────────────────────┤");
    println!("│ Plains (< 0.3):     {:>5.1}%                  │", goals.plains_target);
    println!("│ Hills (0.3-0.6):    {:>5.1}%                  │", goals.hills_target);
    println!("│ Mountains (> 0.6):  {:>5.1}%                  │", goals.mountains_target);
    println!("│ Land Coverage:      {:>5.1}%                  │", goals.land_target);
    println!("│ Tolerance:          ±{:.1}%                    │", goals.tolerance);
    println!("└─────────────────────────────────────────────┘");
    
    // Step 1: Quick grid search to find good starting point
    println!("\n[Phase 1] Quick Grid Search...");
    let _starting_params = quick_search(&goals);
    
    // Step 2: Fine-tune with hill climbing
    println!("\n[Phase 2] Hill-Climbing Optimization...");
    let result = optimize_parameters(&goals, 100, 3);
    
    // Step 3: Print recommendations
    print_recommendations(&result);
    
    // Step 4: Generate a sample map with the optimized parameters to verify
    println!("\n[Phase 3] Verification - Generating sample map...");
    let mut config = GenerationConfig::high_quality(42);
    result.best_params.apply_to_config(&mut config);
    
    let mesh = generate_map_with_config(&config);
    let stats = generate_map_stats(&mesh);
    
    println!("\nVerification Map Statistics:");
    print_map_report(&stats);
    
    // Check if we hit the targets
    let plains_ok = (stats.elevation_stats.plains_percentage - goals.plains_target).abs() < goals.tolerance;
    let hills_ok = (stats.elevation_stats.hills_percentage - goals.hills_target).abs() < goals.tolerance;
    let mountains_ok = (stats.elevation_stats.mountains_percentage - goals.mountains_target).abs() < goals.tolerance;
    
    println!("Target Achievement:");
    println!("  Plains:    {} (target: {:.0}%, actual: {:.1}%)", 
        if plains_ok { "✓" } else { "✗" }, goals.plains_target, stats.elevation_stats.plains_percentage);
    println!("  Hills:     {} (target: {:.0}%, actual: {:.1}%)", 
        if hills_ok { "✓" } else { "✗" }, goals.hills_target, stats.elevation_stats.hills_percentage);
    println!("  Mountains: {} (target: {:.0}%, actual: {:.1}%)", 
        if mountains_ok { "✓" } else { "✗" }, goals.mountains_target, stats.elevation_stats.mountains_percentage);
}
