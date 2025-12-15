//! Elevation calculation for corners and centers.
//!
//! Elevation is computed as "distance from coast" using BFS,
//! then redistributed using a curve to create more flat lowlands.

use std::collections::VecDeque;
use crate::mesh::{DualMesh, NONE};

/// Configuration for elevation generation.
#[derive(Debug, Clone)]
pub struct ElevationConfig {
    /// Redistribution exponent (higher = more lowlands).
    /// Uses formula: y = 1 - (1 - x)^exponent
    pub redistribution_exponent: f64,
    /// Scale factor for lake elevation (lakes are slightly lower than surroundings).
    pub lake_elevation_factor: f64,
}

impl Default for ElevationConfig {
    fn default() -> Self {
        Self {
            // Higher exponent = more lowlands, fewer peaks
            // 2.0 was causing too much high elevation
            redistribution_exponent: 3.0,
            lake_elevation_factor: 0.2,
        }
    }
}

/// Calculate elevation for all corners using BFS from coast.
pub fn calculate_corner_elevation(mesh: &mut DualMesh, config: &ElevationConfig) {
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

    // Normalize and redistribute elevation
    if max_elevation > 0.0 {
        for corner in mesh.corners.iter_mut() {
            if corner.ocean {
                corner.elevation = 0.0;
            } else if corner.water {
                // Lake corners - set to low elevation
                corner.elevation = config.lake_elevation_factor;
            } else if corner.elevation.is_finite() {
                // Normalize to [0, 1]
                let normalized = corner.elevation / max_elevation;
                // Redistribute: y = 1 - (1 - x)^exp
                // This creates more flat lowlands and steeper mountains
                corner.elevation = 1.0 - (1.0 - normalized).powf(config.redistribution_exponent);
            } else {
                // Unreachable corners (shouldn't happen)
                corner.elevation = 1.0;
            }
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
