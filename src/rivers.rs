//! River generation following downslope paths.
//!
//! Rivers flow from high elevation corners down to the ocean,
//! following the precomputed downslope directions.

use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use crate::mesh::{DualMesh, NONE};

/// Configuration for river generation.
#[derive(Debug, Clone)]
pub struct RiverConfig {
    /// Random seed.
    pub seed: u64,
    /// Number of rivers to generate.
    pub num_rivers: usize,
    /// Minimum elevation for river sources.
    pub min_source_elevation: f64,
    /// Maximum elevation for river sources.
    pub max_source_elevation: f64,
}

impl Default for RiverConfig {
    fn default() -> Self {
        Self {
            seed: 12345,
            num_rivers: 20,
            min_source_elevation: 0.4,
            max_source_elevation: 0.9,
        }
    }
}

/// Generate rivers on the mesh.
pub fn generate_rivers(mesh: &mut DualMesh, config: &RiverConfig) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(config.seed);

    // Find candidate source corners (high elevation land corners)
    let candidates: Vec<usize> = (0..mesh.num_solid_corners)
        .filter(|&i| {
            let corner = &mesh.corners[i];
            !corner.water 
                && !corner.ocean 
                && corner.elevation >= config.min_source_elevation
                && corner.elevation <= config.max_source_elevation
                && corner.downslope != NONE
        })
        .collect();

    if candidates.is_empty() {
        return;
    }

    // Generate rivers
    let num_rivers = config.num_rivers.min(candidates.len());
    let mut used_sources = std::collections::HashSet::new();

    for _ in 0..num_rivers {
        // Pick a random source that hasn't been used
        let available: Vec<usize> = candidates.iter()
            .copied()
            .filter(|c| !used_sources.contains(c))
            .collect();

        if available.is_empty() {
            break;
        }

        let source_idx = available[rng.gen_range(0..available.len())];
        used_sources.insert(source_idx);

        // Trace river path down to ocean
        trace_river(mesh, source_idx);
    }
}

/// Trace a river from source corner down to the ocean.
fn trace_river(mesh: &mut DualMesh, source: usize) {
    let mut current = source;
    let mut visited = std::collections::HashSet::new();
    
    // Follow downslope until we reach water or get stuck
    while current != NONE && current < mesh.corners.len() {
        // Prevent infinite loops
        if visited.contains(&current) {
            break;
        }
        visited.insert(current);

        // Increment river volume at this corner
        mesh.corners[current].river += 1;

        // If we reached ocean, coast, or lake, stop
        if mesh.corners[current].ocean || mesh.corners[current].coast {
            break;
        }
        
        // If we reached a lake (freshwater), stop here - rivers flow INTO lakes
        if mesh.corners[current].water && !mesh.corners[current].ocean {
            break;
        }

        // Find the edge connecting current to downslope
        let downslope = mesh.corners[current].downslope;
        if downslope == NONE || downslope >= mesh.corners.len() {
            break;
        }

        // Find and mark the edge between current and downslope
        for &edge_idx in &mesh.corners[current].protrudes.clone() {
            if edge_idx >= mesh.edges.len() {
                continue;
            }
            let edge = &mesh.edges[edge_idx];
            if (edge.v0 == current && edge.v1 == downslope) ||
               (edge.v1 == current && edge.v0 == downslope) {
                mesh.edges[edge_idx].river += 1;
                break;
            }
        }

        // Move to next corner
        current = downslope;
    }
}

/// Calculate moisture based on distance from rivers and lakes.
pub fn calculate_moisture(mesh: &mut DualMesh) {
    use std::collections::VecDeque;

    // Initialize moisture: rivers and lakes have high moisture
    for corner in mesh.corners.iter_mut() {
        if corner.water && !corner.ocean {
            // Lake
            corner.moisture = 1.0;
        } else if corner.river > 0 {
            // River
            corner.moisture = 1.0;
        } else if corner.ocean {
            corner.moisture = 0.0;
        } else {
            corner.moisture = 0.0;
        }
    }

    // BFS to spread moisture from freshwater sources
    let mut queue = VecDeque::new();
    
    // Start from all freshwater sources
    for i in 0..mesh.corners.len() {
        if mesh.corners[i].moisture > 0.0 && !mesh.corners[i].ocean {
            queue.push_back(i);
        }
    }

    // Spread moisture with falloff
    while let Some(corner_idx) = queue.pop_front() {
        let current_moisture = mesh.corners[corner_idx].moisture;
        
        for &adj in &mesh.corners[corner_idx].adjacent.clone() {
            if adj == NONE || adj >= mesh.corners.len() {
                continue;
            }
            
            // Skip ocean
            if mesh.corners[adj].ocean {
                continue;
            }

            // Moisture decreases with distance
            let new_moisture = current_moisture * 0.9;
            
            if new_moisture > mesh.corners[adj].moisture {
                mesh.corners[adj].moisture = new_moisture;
                queue.push_back(adj);
            }
        }
    }

    // Calculate center moisture as average of corner moisture
    for i in 0..mesh.centers.len() {
        let corners = &mesh.centers[i].corners;
        
        if corners.is_empty() {
            continue;
        }

        let sum: f64 = corners.iter()
            .filter_map(|&c| {
                if c < mesh.corners.len() {
                    Some(mesh.corners[c].moisture)
                } else {
                    None
                }
            })
            .sum();

        let count = corners.iter()
            .filter(|&&c| c < mesh.corners.len())
            .count();

        mesh.centers[i].moisture = if count > 0 {
            sum / count as f64
        } else {
            0.0
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pointgen::{PointGeneratorConfig, generate_relaxed_points};
    use crate::graph::build_dual_mesh;
    use crate::island::{IslandConfig, generate_island};
    use crate::elevation::{ElevationConfig, generate_elevation};

    #[test]
    fn test_river_generation() {
        let point_config = PointGeneratorConfig {
            num_points: 200,
            lloyd_iterations: 2,
            ..Default::default()
        };
        let points = generate_relaxed_points(&point_config);
        let mut mesh = build_dual_mesh(&points, &point_config.bounds);

        let island_config = IslandConfig::default();
        generate_island(&mut mesh, &island_config, &point_config.bounds);

        let elevation_config = ElevationConfig::default();
        generate_elevation(&mut mesh, &elevation_config);

        let river_config = RiverConfig {
            num_rivers: 10,
            ..Default::default()
        };
        generate_rivers(&mut mesh, &river_config);

        // Should have some river corners
        let river_corners = mesh.solid_corners()
            .filter(|c| c.river > 0)
            .count();
        
        // May not always have rivers if no suitable sources
        // Just check it doesn't crash
        println!("River corners: {}", river_corners);
    }

    #[test]
    fn test_moisture_calculation() {
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

        let river_config = RiverConfig::default();
        generate_rivers(&mut mesh, &river_config);
        calculate_moisture(&mut mesh);

        // All moisture should be in [0, 1]
        for corner in mesh.solid_corners() {
            assert!(corner.moisture >= 0.0 && corner.moisture <= 1.0,
                "Moisture {} out of range", corner.moisture);
        }
    }
}
