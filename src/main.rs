//! Polygonal Map Generator - Main Application
//!
//! Interactive visualization of the procedural map generation.

use macroquad::prelude::*;
use mapgen::{MapConfig, generate_map};
use mapgen::geometry::BoundingBox;
use mapgen::island::IslandShape;
use mapgen::visualizer::Visualizer;
use mapgen::mesh::DualMesh;

fn window_conf() -> Conf {
    Conf {
        window_title: "Polygonal Map Generator".to_string(),
        window_width: 1280,
        window_height: 800,
        window_resizable: true,
        ..Default::default()
    }
}

struct AppState {
    mesh: DualMesh,
    visualizer: Visualizer,
    config: MapConfig,
    generation_count: u64,
}

impl AppState {
    fn new() -> Self {
        let config = MapConfig {
            bounds: BoundingBox::new(0.0, 0.0, 1000.0, 1000.0),
            num_points: 1000,
            seed: 12345,
            lloyd_iterations: 2,
            island_shape: IslandShape::Radial,
            island_factor: 0.7,
            num_rivers: 30,
        };

        let mesh = generate_map(&config);
        let visualizer = Visualizer::new(config.bounds);

        Self {
            mesh,
            visualizer,
            config,
            generation_count: 0,
        }
    }

    fn regenerate(&mut self) {
        self.generation_count += 1;
        
        // Change seed for variety
        self.config.seed = 12345 + self.generation_count * 1000;
        
        // Cycle through island shapes
        self.config.island_shape = match self.generation_count % 4 {
            0 => IslandShape::Radial,
            1 => IslandShape::Blob,
            2 => IslandShape::Square,
            _ => IslandShape::Noise,
        };

        println!("Regenerating map #{} with {:?} shape...", 
                 self.generation_count, self.config.island_shape);
        
        self.mesh = generate_map(&self.config);
        
        println!("Generated {} regions, {} corners, {} edges",
                 self.mesh.num_solid_centers,
                 self.mesh.num_solid_corners,
                 self.mesh.num_solid_edges);
    }

    fn update(&mut self) {
        // Handle input
        let regenerate = self.visualizer.handle_input();
        
        if regenerate {
            self.regenerate();
        }

        // Handle click selection
        if is_mouse_button_pressed(MouseButton::Left) {
            let mouse_pos = Vec2::new(mouse_position().0, mouse_position().1);
            // Don't select if clicking on UI area
            if mouse_pos.x > 260.0 {
                self.visualizer.selected_center = self.visualizer.find_center_at(&self.mesh, mouse_pos);
            }
        }
    }

    fn draw(&self) {
        self.visualizer.draw(&self.mesh);
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    println!("=== Polygonal Map Generator ===");
    println!("Based on Red Blob Games algorithm");
    println!();
    println!("Controls:");
    println!("  WASD/Arrows - Pan");
    println!("  Mouse Wheel - Zoom");
    println!("  Click - Select region");
    println!("  R - Regenerate map");
    println!("  V - Toggle Voronoi edges");
    println!("  B - Toggle Delaunay edges");
    println!("  F - Toggle rivers");
    println!("  C - Toggle centers");
    println!("  M - Cycle display modes");
    println!("  Home - Reset view");
    println!();

    let mut state = AppState::new();
    
    println!("Initial map generated:");
    println!("  {} regions", state.mesh.num_solid_centers);
    println!("  {} corners", state.mesh.num_solid_corners);
    println!("  {} edges", state.mesh.num_solid_edges);

    loop {
        state.update();
        state.draw();
        next_frame().await
    }
}
