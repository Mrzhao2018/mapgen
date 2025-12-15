//! Macroquad-based visualization for the map generator.
//!
//! Provides interactive debug view with:
//! - Voronoi/Delaunay toggle (V/D keys)
//! - Zoom (mouse wheel) and pan (WASD/arrow keys)
//! - Region selection (click)
//! - Map regeneration (R key)

use macroquad::prelude::*;
use crate::geometry::{Point, BoundingBox};
use crate::mesh::{DualMesh, NONE};
use crate::biomes::biome_color;

/// Display mode for the visualization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DisplayMode {
    /// Show biome colors.
    Biomes,
    /// Show elevation as grayscale.
    Elevation,
    /// Show moisture as blue gradient.
    Moisture,
    /// Show polygon outlines only.
    Wireframe,
}

/// Visualization state.
pub struct Visualizer {
    /// Camera offset (pan).
    pub camera_offset: Vec2,
    /// Zoom level.
    pub zoom: f32,
    /// Show Voronoi edges.
    pub show_voronoi: bool,
    /// Show Delaunay edges.
    pub show_delaunay: bool,
    /// Show rivers.
    pub show_rivers: bool,
    /// Show region centers.
    pub show_centers: bool,
    /// Currently selected center index.
    pub selected_center: Option<usize>,
    /// Display mode.
    pub display_mode: DisplayMode,
    /// Map bounds.
    pub bounds: BoundingBox,
}

impl Default for Visualizer {
    fn default() -> Self {
        Self {
            camera_offset: Vec2::ZERO,
            zoom: 1.0,
            show_voronoi: true,
            show_delaunay: false,
            show_rivers: true,
            show_centers: false,
            selected_center: None,
            display_mode: DisplayMode::Biomes,
            bounds: BoundingBox::new(0.0, 0.0, 1000.0, 1000.0),
        }
    }
}

impl Visualizer {
    pub fn new(bounds: BoundingBox) -> Self {
        // Calculate initial zoom to fit the map
        let screen_w = screen_width();
        let screen_h = screen_height();
        let map_w = bounds.width() as f32;
        let map_h = bounds.height() as f32;
        let zoom = (screen_w / map_w).min(screen_h / map_h) * 0.9;

        // Center the map
        let offset_x = (screen_w - map_w * zoom) / 2.0 - bounds.min.x as f32 * zoom;
        let offset_y = (screen_h - map_h * zoom) / 2.0 - bounds.min.y as f32 * zoom;

        Self {
            camera_offset: Vec2::new(offset_x, offset_y),
            zoom,
            bounds,
            ..Default::default()
        }
    }

    /// Transform a map point to screen coordinates.
    pub fn map_to_screen(&self, p: Point) -> Vec2 {
        Vec2::new(
            p.x as f32 * self.zoom + self.camera_offset.x,
            p.y as f32 * self.zoom + self.camera_offset.y,
        )
    }

    /// Transform screen coordinates to map point.
    pub fn screen_to_map(&self, screen: Vec2) -> Point {
        Point::new(
            ((screen.x - self.camera_offset.x) / self.zoom) as f64,
            ((screen.y - self.camera_offset.y) / self.zoom) as f64,
        )
    }

    /// Handle input for camera control.
    pub fn handle_input(&mut self) -> bool {
        let mut regenerate = false;

        // Zoom with mouse wheel
        let (_, wheel_y) = mouse_wheel();
        if wheel_y != 0.0 {
            let mouse_pos = Vec2::new(mouse_position().0, mouse_position().1);
            let old_map_pos = self.screen_to_map(mouse_pos);
            
            // Reduced sensitivity: 0.05 instead of 0.1
            self.zoom *= 1.0 + wheel_y * 0.05;
            self.zoom = self.zoom.clamp(0.1, 10.0);
            
            // Adjust offset to zoom toward mouse position
            let new_screen_pos = self.map_to_screen(old_map_pos);
            self.camera_offset += mouse_pos - new_screen_pos;
        }

        // Pan with WASD or arrow keys
        let pan_speed = 10.0 / self.zoom;
        if is_key_down(KeyCode::W) || is_key_down(KeyCode::Up) {
            self.camera_offset.y += pan_speed;
        }
        if is_key_down(KeyCode::S) || is_key_down(KeyCode::Down) {
            self.camera_offset.y -= pan_speed;
        }
        if is_key_down(KeyCode::A) || is_key_down(KeyCode::Left) {
            self.camera_offset.x += pan_speed;
        }
        if is_key_down(KeyCode::D) || is_key_down(KeyCode::Right) {
            self.camera_offset.x -= pan_speed;
        }

        // Toggle displays
        if is_key_pressed(KeyCode::V) {
            self.show_voronoi = !self.show_voronoi;
        }
        if is_key_pressed(KeyCode::B) {
            self.show_delaunay = !self.show_delaunay;
        }
        if is_key_pressed(KeyCode::C) {
            self.show_centers = !self.show_centers;
        }
        if is_key_pressed(KeyCode::F) {
            self.show_rivers = !self.show_rivers;
        }

        // Cycle display modes with M
        if is_key_pressed(KeyCode::M) {
            self.display_mode = match self.display_mode {
                DisplayMode::Biomes => DisplayMode::Elevation,
                DisplayMode::Elevation => DisplayMode::Moisture,
                DisplayMode::Moisture => DisplayMode::Wireframe,
                DisplayMode::Wireframe => DisplayMode::Biomes,
            };
        }

        // Regenerate with R
        if is_key_pressed(KeyCode::R) {
            regenerate = true;
        }

        // Reset view with Home
        if is_key_pressed(KeyCode::Home) {
            *self = Self::new(self.bounds);
        }

        regenerate
    }

    /// Find which center contains the given screen position.
    pub fn find_center_at(&self, mesh: &DualMesh, screen_pos: Vec2) -> Option<usize> {
        let map_pos = self.screen_to_map(screen_pos);

        // Simple brute force search
        let mut best_center = None;
        let mut best_distance = f64::INFINITY;

        for (i, center) in mesh.centers.iter().enumerate() {
            if i >= mesh.num_solid_centers {
                continue; // Skip ghost centers
            }
            
            let dist = center.position.distance_sq(&map_pos);
            if dist < best_distance {
                best_distance = dist;
                best_center = Some(i);
            }
        }

        best_center
    }

    /// Draw the complete map.
    pub fn draw(&self, mesh: &DualMesh) {
        // Clear background
        clear_background(Color::from_rgba(20, 30, 40, 255));

        // Draw polygons
        for i in 0..mesh.num_solid_centers {
            self.draw_polygon(mesh, i);
        }

        // Draw Voronoi edges
        if self.show_voronoi {
            self.draw_voronoi_edges(mesh);
        }

        // Draw Delaunay edges
        if self.show_delaunay {
            self.draw_delaunay_edges(mesh);
        }

        // Draw rivers
        if self.show_rivers {
            self.draw_rivers(mesh);
        }

        // Draw centers
        if self.show_centers {
            self.draw_centers(mesh);
        }

        // Draw selected region highlight
        if let Some(idx) = self.selected_center {
            self.draw_selected(mesh, idx);
        }

        // Draw UI
        self.draw_ui(mesh);
    }

    /// Draw a single Voronoi polygon.
    fn draw_polygon(&self, mesh: &DualMesh, center_idx: usize) {
        let center = &mesh.centers[center_idx];
        // Use clipped vertices to avoid drawing outside bounds
        let vertices = mesh.get_polygon_vertices_clipped(center_idx, Some(&self.bounds));

        if vertices.len() < 3 {
            return;
        }

        // Get color based on display mode
        let color = match self.display_mode {
            DisplayMode::Biomes => {
                let (r, g, b) = biome_color(center.biome);
                Color::from_rgba(r, g, b, 255)
            }
            DisplayMode::Elevation => {
                let e = (center.elevation * 255.0) as u8;
                Color::from_rgba(e, e, e, 255)
            }
            DisplayMode::Moisture => {
                let m = (center.moisture * 255.0) as u8;
                Color::from_rgba(0, m / 2, m, 255)
            }
            DisplayMode::Wireframe => {
                Color::from_rgba(30, 40, 50, 255)
            }
        };

        // Convert to screen coordinates
        let screen_vertices: Vec<Vec2> = vertices.iter()
            .map(|v| self.map_to_screen(*v))
            .collect();

        // Draw filled polygon using triangle fan
        if screen_vertices.len() >= 3 {
            let center_screen = self.map_to_screen(center.position);
            
            for i in 0..screen_vertices.len() {
                let v1 = screen_vertices[i];
                let v2 = screen_vertices[(i + 1) % screen_vertices.len()];
                draw_triangle(center_screen, v1, v2, color);
            }
        }
    }

    /// Draw Voronoi edges.
    fn draw_voronoi_edges(&self, mesh: &DualMesh) {
        let color = Color::from_rgba(60, 60, 60, 255);
        let thickness = 1.0;

        for edge in mesh.solid_edges() {
            if edge.v0 == NONE || edge.v1 == NONE {
                continue;
            }
            if edge.v0 >= mesh.corners.len() || edge.v1 >= mesh.corners.len() {
                continue;
            }

            // Clamp corner positions to bounds
            let c0 = &mesh.corners[edge.v0].position;
            let c1 = &mesh.corners[edge.v1].position;
            let pos0 = Point::new(
                c0.x.clamp(self.bounds.min.x, self.bounds.max.x),
                c0.y.clamp(self.bounds.min.y, self.bounds.max.y),
            );
            let pos1 = Point::new(
                c1.x.clamp(self.bounds.min.x, self.bounds.max.x),
                c1.y.clamp(self.bounds.min.y, self.bounds.max.y),
            );

            let p0 = self.map_to_screen(pos0);
            let p1 = self.map_to_screen(pos1);
            draw_line(p0.x, p0.y, p1.x, p1.y, thickness, color);
        }
    }

    /// Draw Delaunay edges.
    fn draw_delaunay_edges(&self, mesh: &DualMesh) {
        let color = Color::from_rgba(100, 50, 50, 255);
        let thickness = 1.0;

        for edge in mesh.solid_edges() {
            if edge.d0 == NONE || edge.d1 == NONE {
                continue;
            }
            if edge.d0 >= mesh.centers.len() || edge.d1 >= mesh.centers.len() {
                continue;
            }

            let p0 = self.map_to_screen(mesh.centers[edge.d0].position);
            let p1 = self.map_to_screen(mesh.centers[edge.d1].position);
            draw_line(p0.x, p0.y, p1.x, p1.y, thickness, color);
        }
    }

    /// Draw rivers.
    fn draw_rivers(&self, mesh: &DualMesh) {
        // Draw river edges - iterate all edges, not just solid
        for edge in &mesh.edges {
            if edge.river == 0 {
                continue;
            }
            if edge.v0 == NONE || edge.v1 == NONE {
                continue;
            }
            if edge.v0 >= mesh.corners.len() || edge.v1 >= mesh.corners.len() {
                continue;
            }

            // Clamp to bounds
            let c0 = &mesh.corners[edge.v0].position;
            let c1 = &mesh.corners[edge.v1].position;
            let pos0 = Point::new(
                c0.x.clamp(self.bounds.min.x, self.bounds.max.x),
                c0.y.clamp(self.bounds.min.y, self.bounds.max.y),
            );
            let pos1 = Point::new(
                c1.x.clamp(self.bounds.min.x, self.bounds.max.x),
                c1.y.clamp(self.bounds.min.y, self.bounds.max.y),
            );

            let p0 = self.map_to_screen(pos0);
            let p1 = self.map_to_screen(pos1);
            
            // Thickness based on river volume - minimum 2px, scale with zoom
            let base_thickness = 1.5 + (edge.river as f32).sqrt() * 1.5;
            let thickness = (base_thickness * self.zoom).max(1.5);
            // Brighter blue color for better visibility
            let color = Color::from_rgba(100, 149, 237, 255); // Cornflower blue
            draw_line(p0.x, p0.y, p1.x, p1.y, thickness, color);
        }
    }

    /// Draw center points.
    fn draw_centers(&self, mesh: &DualMesh) {
        for i in 0..mesh.num_solid_centers {
            let center = &mesh.centers[i];
            let pos = self.map_to_screen(center.position);
            let radius = 3.0;
            let color = if center.water {
                Color::from_rgba(51, 102, 153, 255)
            } else {
                Color::from_rgba(200, 100, 100, 255)
            };
            draw_circle(pos.x, pos.y, radius, color);
        }
    }

    /// Draw selected region highlight.
    fn draw_selected(&self, mesh: &DualMesh, center_idx: usize) {
        if center_idx >= mesh.centers.len() {
            return;
        }

        // Use clipped vertices
        let vertices = mesh.get_polygon_vertices_clipped(center_idx, Some(&self.bounds));
        if vertices.len() < 3 {
            return;
        }

        let screen_vertices: Vec<Vec2> = vertices.iter()
            .map(|v| self.map_to_screen(*v))
            .collect();

        // Draw highlight border
        let color = Color::from_rgba(255, 255, 0, 255);
        let thickness = 3.0;

        for i in 0..screen_vertices.len() {
            let v1 = screen_vertices[i];
            let v2 = screen_vertices[(i + 1) % screen_vertices.len()];
            draw_line(v1.x, v1.y, v2.x, v2.y, thickness, color);
        }

        // Draw center point
        let center_pos = self.map_to_screen(mesh.centers[center_idx].position);
        draw_circle(center_pos.x, center_pos.y, 5.0, color);
    }

    /// Draw UI overlay.
    fn draw_ui(&self, mesh: &DualMesh) {
        let mut y = 10.0;
        let line_height = 18.0;
        let font_size = 16.0;

        // Title
        draw_text("Polygonal Map Generator", 10.0, y, font_size + 4.0, WHITE);
        y += line_height + 5.0;

        // Controls
        draw_text("Controls:", 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text("  WASD/Arrows: Pan", 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text("  Mouse Wheel: Zoom", 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text("  Click: Select region", 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text("  R: Regenerate", 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text("  Home: Reset view", 10.0, y, font_size, GRAY);
        y += line_height + 5.0;

        // Toggles
        draw_text("Toggles:", 10.0, y, font_size, GRAY);
        y += line_height;
        let voronoi_text = format!("  V: Voronoi [{}]", if self.show_voronoi { "ON" } else { "OFF" });
        draw_text(&voronoi_text, 10.0, y, font_size, if self.show_voronoi { GREEN } else { GRAY });
        y += line_height;
        let delaunay_text = format!("  B: Delaunay [{}]", if self.show_delaunay { "ON" } else { "OFF" });
        draw_text(&delaunay_text, 10.0, y, font_size, if self.show_delaunay { GREEN } else { GRAY });
        y += line_height;
        let rivers_text = format!("  F: Rivers [{}]", if self.show_rivers { "ON" } else { "OFF" });
        draw_text(&rivers_text, 10.0, y, font_size, if self.show_rivers { GREEN } else { GRAY });
        y += line_height;
        let centers_text = format!("  C: Centers [{}]", if self.show_centers { "ON" } else { "OFF" });
        draw_text(&centers_text, 10.0, y, font_size, if self.show_centers { GREEN } else { GRAY });
        y += line_height;
        let mode_name = match self.display_mode {
            DisplayMode::Biomes => "Biomes",
            DisplayMode::Elevation => "Elevation",
            DisplayMode::Moisture => "Moisture",
            DisplayMode::Wireframe => "Wireframe",
        };
        draw_text(&format!("  M: Mode [{}]", mode_name), 10.0, y, font_size, YELLOW);
        y += line_height + 5.0;

        // Stats
        draw_text("Stats:", 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text(&format!("  Regions: {}", mesh.num_solid_centers), 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text(&format!("  Corners: {}", mesh.num_solid_corners), 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text(&format!("  Edges: {}", mesh.num_solid_edges), 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text(&format!("  Zoom: {:.1}x", self.zoom), 10.0, y, font_size, GRAY);
        y += line_height;
        draw_text(&format!("  FPS: {}", get_fps()), 10.0, y, font_size, GRAY);

        // Selected region info
        if let Some(idx) = self.selected_center {
            if idx < mesh.centers.len() {
                let center = &mesh.centers[idx];
                
                let info_x = screen_width() - 250.0;
                let mut info_y = 10.0;
                
                draw_text("Selected Region:", info_x, info_y, font_size + 2.0, YELLOW);
                info_y += line_height;
                draw_text(&format!("  Index: {}", idx), info_x, info_y, font_size, WHITE);
                info_y += line_height;
                draw_text(&format!("  Biome: {:?}", center.biome), info_x, info_y, font_size, WHITE);
                info_y += line_height;
                draw_text(&format!("  Elevation: {:.3}", center.elevation), info_x, info_y, font_size, WHITE);
                info_y += line_height;
                draw_text(&format!("  Moisture: {:.3}", center.moisture), info_x, info_y, font_size, WHITE);
                info_y += line_height;
                draw_text(&format!("  Water: {}", center.water), info_x, info_y, font_size, WHITE);
                info_y += line_height;
                draw_text(&format!("  Ocean: {}", center.ocean), info_x, info_y, font_size, WHITE);
                info_y += line_height;
                draw_text(&format!("  Coast: {}", center.coast), info_x, info_y, font_size, WHITE);
                info_y += line_height;
                draw_text(&format!("  Neighbors: {}", center.neighbors.len()), info_x, info_y, font_size, WHITE);
                info_y += line_height;
                draw_text(&format!("  Corners: {}", center.corners.len()), info_x, info_y, font_size, WHITE);
            }
        }
    }
}
