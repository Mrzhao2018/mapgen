//! Elevation calculation for corners and centers.
//!
//! Elevation is computed using a hybrid approach:
//! 1. BFS "distance from coast" as base gradient
//! 2. Noise modulation (Fbm/RidgedMulti) for natural variation
//! 3. Pit-filling algorithm to eliminate local depressions
//! 4. Redistribution curve for more flat lowlands

use std::collections::VecDeque;
use noise::{NoiseFn, Fbm, RidgedMulti, Perlin, MultiFractal};
use crate::mesh::{DualMesh, NONE};

/// Configuration for elevation generation.
#[derive(Debug, Clone)]
pub struct ElevationConfig {
    /// Redistribution exponent (higher = more lowlands).
    /// Uses formula: y = 1 - (1 - x)^exponent
    pub redistribution_exponent: f64,
    /// Scale factor for lake elevation (lakes are slightly lower than surroundings).
    pub lake_elevation_factor: f64,
    /// Use hybrid elevation (BFS + noise modulation)
    pub use_hybrid: bool,
    /// Weight for BFS component in hybrid mode [0, 1]
    pub bfs_weight: f64,
    /// Noise frequency for elevation modulation
    pub noise_frequency: f64,
    /// Enable pit-filling algorithm
    pub fill_pits: bool,
    /// Seed for noise generation
    pub seed: u32,
    /// Use ridged multi-fractal instead of Fbm
    pub use_ridged: bool,
}

impl Default for ElevationConfig {
    fn default() -> Self {
        Self {
            // Optimized for ~40% plains, ~35% hills, ~25% mountains
            redistribution_exponent: 0.75,  // Lower = more mountains (0.75 for ~25% mountains)
            lake_elevation_factor: 0.08,
            use_hybrid: true,
            bfs_weight: 0.55,  // Lower = more noise influence = more varied terrain
            noise_frequency: 0.002,
            fill_pits: true,
            seed: 42,
            use_ridged: true,
        }
    }
}

/// Noise-based elevation generator.
pub struct ElevationNoise {
    fbm: Fbm<Perlin>,
    ridged: RidgedMulti<Perlin>,
    use_ridged: bool,
}

impl ElevationNoise {
    pub fn new(config: &ElevationConfig) -> Self {
        let mut fbm = Fbm::<Perlin>::new(config.seed);
        fbm = fbm
            .set_octaves(3)
            .set_frequency(config.noise_frequency)
            .set_lacunarity(2.5)
            .set_persistence(0.4);

        let mut ridged = RidgedMulti::<Perlin>::new(config.seed);
        ridged = ridged
            .set_octaves(3)
            .set_frequency(config.noise_frequency * 0.4)  // Very low freq for broad ridges
            .set_lacunarity(2.5)
            .set_persistence(0.65);  // Higher persistence for sharper ridges

        Self {
            fbm,
            ridged,
            use_ridged: config.use_ridged,
        }
    }

    /// Sample elevation noise at a point, returning value in [0, 1]
    pub fn sample(&self, x: f64, y: f64) -> f64 {
        let raw = if self.use_ridged {
            // RidgedMulti produces mountain-like ridges
            // Output is typically in [0, 1] range
            self.ridged.get([x, y])
        } else {
            // Fbm output is typically in [-1, 1] range
            (self.fbm.get([x, y]) + 1.0) * 0.5
        };
        raw.clamp(0.0, 1.0)
    }
}

/// Calculate elevation for all corners using BFS from coast.
pub fn calculate_corner_elevation(mesh: &mut DualMesh, config: &ElevationConfig) {
    // Create noise generator if using hybrid mode
    let noise = if config.use_hybrid {
        Some(ElevationNoise::new(config))
    } else {
        None
    };

    // Initialize all corners
    for corner in mesh.corners.iter_mut() {
        if corner.ocean {
            corner.elevation = 0.0;
        } else if corner.water {
            // Lakes - will be set later based on surrounding land
            corner.elevation = 0.0;
        } else {
            corner.elevation = f64::INFINITY;
        }
    }

    // BFS from coast corners
    let mut queue = VecDeque::new();
    
    // Start from coast corners (adjacent to ocean)
    for i in 0..mesh.corners.len() {
        if mesh.corners[i].coast && !mesh.corners[i].ocean {
            mesh.corners[i].elevation = 0.0;
            queue.push_back(i);
        }
    }

    // Also start from land corners adjacent to ocean
    for i in 0..mesh.corners.len() {
        if !mesh.corners[i].water && !mesh.corners[i].ocean {
            let has_ocean_neighbor = mesh.corners[i].adjacent.iter()
                .any(|&adj| adj != NONE && adj < mesh.corners.len() && mesh.corners[adj].ocean);
            
            if has_ocean_neighbor {
                mesh.corners[i].elevation = 0.0;
                if !queue.contains(&i) {
                    queue.push_back(i);
                }
            }
        }
    }

    // BFS to propagate elevation inland
    while let Some(corner_idx) = queue.pop_front() {
        let current_elevation = mesh.corners[corner_idx].elevation;
        
        for &adj in &mesh.corners[corner_idx].adjacent.clone() {
            if adj == NONE || adj >= mesh.corners.len() {
                continue;
            }
            
            // Skip ocean corners
            if mesh.corners[adj].ocean {
                continue;
            }

            let new_elevation = current_elevation + 1.0;
            
            if new_elevation < mesh.corners[adj].elevation {
                mesh.corners[adj].elevation = new_elevation;
                queue.push_back(adj);
            }
        }
    }

    // Find maximum elevation for normalization
    let max_elevation = mesh.corners.iter()
        .filter(|c| !c.water && c.elevation.is_finite())
        .map(|c| c.elevation)
        .fold(0.0f64, |a, b| a.max(b));

    // Normalize and apply hybrid elevation
    if max_elevation > 0.0 {
        for i in 0..mesh.corners.len() {
            let corner = &mesh.corners[i];
            
            if corner.ocean {
                mesh.corners[i].elevation = 0.0;
            } else if corner.water {
                // Lake corners - set to low elevation
                mesh.corners[i].elevation = config.lake_elevation_factor;
            } else if corner.elevation.is_finite() {
                // Normalize BFS elevation to [0, 1]
                let bfs_elevation = corner.elevation / max_elevation;
                
                // Apply hybrid mode: blend BFS with noise
                let final_elevation = if let Some(ref noise_gen) = noise {
                    let pos = corner.position;
                    let noise_value = noise_gen.sample(pos.x, pos.y);
                    
                    // BFS provides the main gradient, noise only affects high areas
                    // This keeps plains flat and only adds variation to mountains
                    let bfs_component = bfs_elevation;
                    let noise_weight = (1.0 - config.bfs_weight) * bfs_elevation.powf(2.0);
                    
                    // Noise only kicks in at higher elevations
                    let blended = bfs_component + noise_value * noise_weight * 0.3;
                    blended.clamp(0.0, 1.0)
                } else {
                    bfs_elevation
                };
                
                // Redistribute using power function: y = x^exp
                // Higher exponent = more flat lowlands (pushes values toward 0)
                // exp > 1 creates more plains, exp < 1 creates more mountains
                mesh.corners[i].elevation = final_elevation.powf(config.redistribution_exponent);
            } else {
                // Unreachable corners (shouldn't happen)
                mesh.corners[i].elevation = 1.0;
            }
        }
    }

    // Apply pit-filling if enabled
    if config.fill_pits {
        fill_pits(mesh);
    }
}

/// Fill pits (local depressions) in the elevation field.
/// Uses the Planchon-Darboux algorithm for pit removal.
/// This ensures water can always drain to the ocean.
pub fn fill_pits(mesh: &mut DualMesh) {
    const EPSILON: f64 = 0.0001;
    const MAX_ITERATIONS: usize = 100;
    
    // We need to fill any local minimum that isn't ocean
    // A pit is a corner that is lower than all its neighbors but not ocean
    
    for _ in 0..MAX_ITERATIONS {
        let mut changed = false;
        
        for i in 0..mesh.corners.len() {
            // Skip ocean corners - they are valid sinks
            if mesh.corners[i].ocean || mesh.corners[i].water {
                continue;
            }
            
            // Check if this is a local minimum (pit)
            let current_elevation = mesh.corners[i].elevation;
            let mut min_neighbor_elevation = f64::INFINITY;
            let mut has_lower_neighbor = false;
            let mut has_ocean_neighbor = false;
            
            let adjacent = mesh.corners[i].adjacent.clone();
            for &adj in &adjacent {
                if adj == NONE || adj >= mesh.corners.len() {
                    continue;
                }
                
                let adj_elevation = mesh.corners[adj].elevation;
                
                if mesh.corners[adj].ocean {
                    has_ocean_neighbor = true;
                }
                
                if adj_elevation < current_elevation {
                    has_lower_neighbor = true;
                }
                
                if adj_elevation < min_neighbor_elevation {
                    min_neighbor_elevation = adj_elevation;
                }
            }
            
            // If no lower neighbor and not draining to ocean, this is a pit
            if !has_lower_neighbor && !has_ocean_neighbor && min_neighbor_elevation.is_finite() {
                // Raise elevation to slightly above lowest neighbor
                // This allows water to flow out
                let new_elevation = min_neighbor_elevation + EPSILON;
                if new_elevation > current_elevation {
                    mesh.corners[i].elevation = new_elevation.min(1.0);
                    changed = true;
                }
            }
        }
        
        if !changed {
            break;
        }
    }
}

/// Calculate elevation for centers based on their corners.
pub fn calculate_center_elevation(mesh: &mut DualMesh) {
    for i in 0..mesh.centers.len() {
        let corners = &mesh.centers[i].corners;
        
        if corners.is_empty() {
            mesh.centers[i].elevation = 0.0;
            continue;
        }

        if mesh.centers[i].ocean {
            mesh.centers[i].elevation = 0.0;
            continue;
        }

        // Average elevation of corners
        let sum: f64 = corners.iter()
            .filter_map(|&c| {
                if c < mesh.corners.len() {
                    Some(mesh.corners[c].elevation)
                } else {
                    None
                }
            })
            .sum();

        let count = corners.iter()
            .filter(|&&c| c < mesh.corners.len())
            .count();

        mesh.centers[i].elevation = if count > 0 {
            sum / count as f64
        } else {
            0.0
        };
    }
}

/// Calculate downslope direction for each corner.
/// Each corner points to its lowest adjacent neighbor.
pub fn calculate_downslope(mesh: &mut DualMesh) {
    for i in 0..mesh.corners.len() {
        let current_elevation = mesh.corners[i].elevation;
        let mut lowest_neighbor = NONE;
        let mut lowest_elevation = current_elevation;

        for &adj in &mesh.corners[i].adjacent.clone() {
            if adj == NONE || adj >= mesh.corners.len() {
                continue;
            }

            let adj_elevation = mesh.corners[adj].elevation;
            if adj_elevation < lowest_elevation {
                lowest_elevation = adj_elevation;
                lowest_neighbor = adj;
            }
        }

        mesh.corners[i].downslope = lowest_neighbor;
    }
}

/// Generate elevation for the entire mesh.
pub fn generate_elevation(mesh: &mut DualMesh, config: &ElevationConfig) {
    calculate_corner_elevation(mesh, config);
    calculate_center_elevation(mesh);
    calculate_downslope(mesh);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pointgen::{PointGeneratorConfig, generate_relaxed_points};
    use crate::graph::build_dual_mesh;
    use crate::island::{IslandConfig, generate_island};

    #[test]
    fn test_elevation_generation() {
        let point_config = PointGeneratorConfig {
            num_points: 100,
            lloyd_iterations: 2,
            ..Default::default()
        };
        let points = generate_relaxed_points(&point_config);
        let mut mesh = build_dual_mesh(&points, &point_config.bounds);

        let island_config = IslandConfig::default();
        generate_island(&mut mesh, &island_config, &point_config.bounds);

        let elevation_config = ElevationConfig::default();
        generate_elevation(&mut mesh, &elevation_config);

        // Ocean corners should have 0 elevation
        for corner in mesh.solid_corners() {
            if corner.ocean {
                assert_eq!(corner.elevation, 0.0, "Ocean corner should have 0 elevation");
            }
        }

        // Land corners should have positive elevation
        let has_positive_elevation = mesh.solid_corners()
            .any(|c| !c.water && c.elevation > 0.0);
        assert!(has_positive_elevation, "Should have some land with positive elevation");

        // All elevations should be in [0, 1]
        for corner in mesh.solid_corners() {
            assert!(corner.elevation >= 0.0 && corner.elevation <= 1.0,
                "Elevation {} out of range", corner.elevation);
        }
    }

    #[test]
    fn test_downslope() {
        let point_config = PointGeneratorConfig {
            num_points: 50,
            lloyd_iterations: 1,
            ..Default::default()
        };
        let points = generate_relaxed_points(&point_config);
        let mut mesh = build_dual_mesh(&points, &point_config.bounds);

        let island_config = IslandConfig::default();
        generate_island(&mut mesh, &island_config, &point_config.bounds);

        let elevation_config = ElevationConfig::default();
        generate_elevation(&mut mesh, &elevation_config);

        // Downslope should point to lower elevation or NONE
        for corner in mesh.solid_corners() {
            if corner.downslope != NONE && corner.downslope < mesh.corners.len() {
                let downslope_elevation = mesh.corners[corner.downslope].elevation;
                assert!(downslope_elevation <= corner.elevation,
                    "Downslope should point to lower or equal elevation");
            }
        }
    }
}
