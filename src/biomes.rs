//! Biome assignment based on elevation and moisture.
//!
//! Uses a Whittaker diagram to map (elevation, moisture) to biome types.

use crate::mesh::{DualMesh, Biome};

/// Assign biomes to all centers based on elevation and moisture.
pub fn assign_biomes(mesh: &mut DualMesh) {
    for i in 0..mesh.centers.len() {
        let center = &mesh.centers[i];
        
        let biome = if center.ocean {
            Biome::Ocean
        } else if center.water {
            // Freshwater lake
            if center.elevation > 0.8 {
                Biome::Ice
            } else if center.elevation < 0.1 {
                Biome::Marsh
            } else {
                Biome::Lake
            }
        } else if center.coast {
            Biome::Beach
        } else {
            // Land - use Whittaker diagram
            get_biome(center.elevation, center.moisture)
        };

        mesh.centers[i].biome = biome;
    }
}

/// Get biome based on Whittaker diagram.
/// Elevation: 0.0 (sea level) to 1.0 (mountain peak)
/// Moisture: 0.0 (dry) to 1.0 (wet)
fn get_biome(elevation: f64, moisture: f64) -> Biome {
    // High elevation biomes (peaks only)
    if elevation > 0.9 {
        if moisture > 0.5 {
            return Biome::Snow;
        } else if moisture > 0.33 {
            return Biome::Tundra;
        } else if moisture > 0.16 {
            return Biome::Bare;
        } else {
            return Biome::Scorched;
        }
    }

    // Medium-high elevation (mountains)
    if elevation > 0.7 {
        if moisture > 0.66 {
            return Biome::Taiga;
        } else if moisture > 0.33 {
            return Biome::Shrubland;
        } else {
            return Biome::TemperateDesert;
        }
    }

    // Medium elevation (hills)
    if elevation > 0.4 {
        if moisture > 0.83 {
            return Biome::TemperateRainForest;
        } else if moisture > 0.5 {
            return Biome::TemperateDeciduousForest;
        } else if moisture > 0.16 {
            return Biome::Grassland;
        } else {
            return Biome::TemperateDesert;
        }
    }

    // Low elevation
    if moisture > 0.66 {
        return Biome::TropicalRainForest;
    } else if moisture > 0.33 {
        return Biome::TropicalSeasonalForest;
    } else if moisture > 0.16 {
        return Biome::Grassland;
    } else {
        return Biome::SubtropicalDesert;
    }
}

/// Get the display color for a biome (RGB).
pub fn biome_color(biome: Biome) -> (u8, u8, u8) {
    match biome {
        Biome::Ocean => (68, 68, 122),
        Biome::Lake => (51, 102, 153),
        Biome::Ice => (220, 220, 255),
        Biome::Marsh => (47, 102, 102),
        Biome::Beach => (194, 178, 128),
        Biome::Snow => (248, 248, 248),
        Biome::Tundra => (187, 187, 170),
        Biome::Bare => (136, 136, 136),
        Biome::Scorched => (85, 85, 85),
        Biome::Taiga => (153, 170, 119),
        Biome::Shrubland => (136, 153, 119),
        Biome::TemperateDesert => (201, 210, 155),
        Biome::TemperateRainForest => (68, 136, 85),
        Biome::TemperateDeciduousForest => (103, 148, 89),
        Biome::Grassland => (136, 170, 85),
        Biome::TropicalRainForest => (51, 119, 85),
        Biome::TropicalSeasonalForest => (85, 153, 68),
        Biome::SubtropicalDesert => (210, 185, 139),
    }
}

/// Get the display name for a biome.
pub fn biome_name(biome: Biome) -> &'static str {
    match biome {
        Biome::Ocean => "Ocean",
        Biome::Lake => "Lake",
        Biome::Ice => "Ice",
        Biome::Marsh => "Marsh",
        Biome::Beach => "Beach",
        Biome::Snow => "Snow",
        Biome::Tundra => "Tundra",
        Biome::Bare => "Bare",
        Biome::Scorched => "Scorched",
        Biome::Taiga => "Taiga",
        Biome::Shrubland => "Shrubland",
        Biome::TemperateDesert => "Temperate Desert",
        Biome::TemperateRainForest => "Temperate Rain Forest",
        Biome::TemperateDeciduousForest => "Temperate Deciduous Forest",
        Biome::Grassland => "Grassland",
        Biome::TropicalRainForest => "Tropical Rain Forest",
        Biome::TropicalSeasonalForest => "Tropical Seasonal Forest",
        Biome::SubtropicalDesert => "Subtropical Desert",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pointgen::{PointGeneratorConfig, generate_relaxed_points};
    use crate::graph::build_dual_mesh;
    use crate::island::{IslandConfig, generate_island};
    use crate::elevation::{ElevationConfig, generate_elevation};
    use crate::rivers::{RiverConfig, generate_rivers, calculate_moisture};

    #[test]
    fn test_biome_assignment() {
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

        let river_config = RiverConfig::default();
        generate_rivers(&mut mesh, &river_config);
        calculate_moisture(&mut mesh);

        assign_biomes(&mut mesh);

        // Should have variety of biomes
        let mut biome_counts = std::collections::HashMap::new();
        for center in mesh.solid_centers() {
            *biome_counts.entry(center.biome).or_insert(0) += 1;
        }

        // Should have at least ocean
        assert!(biome_counts.contains_key(&Biome::Ocean), "Should have ocean biome");
    }

    #[test]
    fn test_biome_colors() {
        // All biomes should have valid colors
        let biomes = [
            Biome::Ocean, Biome::Lake, Biome::Ice, Biome::Marsh, Biome::Beach,
            Biome::Snow, Biome::Tundra, Biome::Bare, Biome::Scorched, Biome::Taiga,
            Biome::Shrubland, Biome::TemperateDesert, Biome::TemperateRainForest,
            Biome::TemperateDeciduousForest, Biome::Grassland, Biome::TropicalRainForest,
            Biome::TropicalSeasonalForest, Biome::SubtropicalDesert,
        ];

        for biome in biomes {
            let (_r, _g, _b) = biome_color(biome);
            let _name = biome_name(biome);
            // Just ensure no panic
        }
    }
}
