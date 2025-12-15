//! Point generation and Lloyd relaxation.

use crate::geometry::{Point, BoundingBox, polygon_centroid};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Configuration for point generation.
#[derive(Debug, Clone)]
pub struct PointGeneratorConfig {
    /// Bounding box for the map.
    pub bounds: BoundingBox,
    /// Number of points to generate.
    pub num_points: usize,
    /// Random seed for reproducibility.
    pub seed: u64,
    /// Number of Lloyd relaxation iterations.
    pub lloyd_iterations: u32,
    /// Over-relaxation factor (1.0 = standard, 1.5-1.8 for faster convergence).
    pub lloyd_omega: f64,
}

impl Default for PointGeneratorConfig {
    fn default() -> Self {
        Self {
            bounds: BoundingBox::new(0.0, 0.0, 1000.0, 1000.0),
            num_points: 1000,
            seed: 12345,
            lloyd_iterations: 2,
            lloyd_omega: 1.0,
        }
    }
}

/// Generate random points within the bounding box.
pub fn generate_random_points(config: &PointGeneratorConfig) -> Vec<Point> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(config.seed);
    let bounds = &config.bounds;

    let mut points = Vec::with_capacity(config.num_points);
    for _ in 0..config.num_points {
        let x = rng.gen_range(bounds.min.x..bounds.max.x);
        let y = rng.gen_range(bounds.min.y..bounds.max.y);
        points.push(Point::new(x, y));
    }

    points
}

/// Perform one iteration of Lloyd relaxation.
/// Returns the maximum displacement of any point.
fn lloyd_iteration(points: &mut [Point], bounds: &BoundingBox, omega: f64) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }

    // Convert to delaunator points
    let delaunator_points: Vec<delaunator::Point> = points
        .iter()
        .map(|p| delaunator::Point { x: p.x, y: p.y })
        .collect();

    // Triangulate
    let triangulation = delaunator::triangulate(&delaunator_points);
    if triangulation.triangles.is_empty() {
        return 0.0;
    }

    // Build a mapping from point index to one of its incoming half-edges
    let mut point_to_halfedge: Vec<usize> = vec![delaunator::EMPTY; points.len()];
    for e in 0..triangulation.triangles.len() {
        let endpoint = triangulation.triangles[next_halfedge(e)];
        if point_to_halfedge[endpoint] == delaunator::EMPTY || triangulation.halfedges[e] == delaunator::EMPTY {
            point_to_halfedge[endpoint] = e;
        }
    }

    // Compute the centroid of each Voronoi cell
    let mut new_positions = vec![Point::ZERO; points.len()];
    let mut max_displacement = 0.0f64;

    for p in 0..points.len() {
        let start = point_to_halfedge[p];
        if start == delaunator::EMPTY {
            new_positions[p] = points[p];
            continue;
        }

        // Collect triangles around this point
        let mut triangles = Vec::new();
        let mut incoming = start;
        loop {
            triangles.push(incoming / 3);
            let outgoing = next_halfedge(incoming);
            incoming = triangulation.halfedges[outgoing];
            if incoming == delaunator::EMPTY || incoming == start {
                break;
            }
        }

        if triangles.is_empty() {
            new_positions[p] = points[p];
            continue;
        }

        // Compute circumcenters (Voronoi vertices) for these triangles
        let mut voronoi_vertices = Vec::with_capacity(triangles.len());
        for &t in &triangles {
            let i0 = triangulation.triangles[3 * t];
            let i1 = triangulation.triangles[3 * t + 1];
            let i2 = triangulation.triangles[3 * t + 2];

            let p0 = &points[i0];
            let p1 = &points[i1];
            let p2 = &points[i2];

            let cc = crate::geometry::circumcenter(p0, p1, p2);
            
            // Clamp to bounds
            let clamped = Point::new(
                cc.x.clamp(bounds.min.x, bounds.max.x),
                cc.y.clamp(bounds.min.y, bounds.max.y),
            );
            voronoi_vertices.push(clamped);
        }

        // Compute centroid of the Voronoi cell
        let centroid = polygon_centroid(&voronoi_vertices);
        
        // Clamp centroid to bounds
        let centroid = Point::new(
            centroid.x.clamp(bounds.min.x, bounds.max.x),
            centroid.y.clamp(bounds.min.y, bounds.max.y),
        );

        // Apply over-relaxation
        let old_pos = points[p];
        let new_pos = Point::new(
            old_pos.x + omega * (centroid.x - old_pos.x),
            old_pos.y + omega * (centroid.y - old_pos.y),
        );

        // Clamp final position
        new_positions[p] = Point::new(
            new_pos.x.clamp(bounds.min.x, bounds.max.x),
            new_pos.y.clamp(bounds.min.y, bounds.max.y),
        );

        let displacement = old_pos.distance(&new_positions[p]);
        max_displacement = max_displacement.max(displacement);
    }

    // Update points
    points.copy_from_slice(&new_positions);
    max_displacement
}

/// Apply Lloyd relaxation to improve point distribution.
pub fn lloyd_relaxation(points: &mut [Point], bounds: &BoundingBox, iterations: u32, omega: f64) {
    for _ in 0..iterations {
        let _displacement = lloyd_iteration(points, bounds, omega);
    }
}

/// Generate points and apply Lloyd relaxation.
pub fn generate_relaxed_points(config: &PointGeneratorConfig) -> Vec<Point> {
    let mut points = generate_random_points(config);
    lloyd_relaxation(&mut points, &config.bounds, config.lloyd_iterations, config.lloyd_omega);
    points
}

// === Halfedge navigation helpers ===

/// Get the next halfedge in a triangle (CCW).
#[inline]
fn next_halfedge(e: usize) -> usize {
    if e % 3 == 2 { e - 2 } else { e + 1 }
}

/// Get the previous halfedge in a triangle (CW).
#[inline]
#[allow(dead_code)]
fn prev_halfedge(e: usize) -> usize {
    if e % 3 == 0 { e + 2 } else { e - 1 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_random_points() {
        let config = PointGeneratorConfig {
            num_points: 100,
            ..Default::default()
        };
        let points = generate_random_points(&config);
        assert_eq!(points.len(), 100);

        // All points should be within bounds
        for p in &points {
            assert!(config.bounds.contains(p));
        }
    }

    #[test]
    fn test_deterministic_generation() {
        let config = PointGeneratorConfig {
            num_points: 50,
            seed: 42,
            ..Default::default()
        };
        let points1 = generate_random_points(&config);
        let points2 = generate_random_points(&config);
        
        // Same seed should produce same points
        assert_eq!(points1.len(), points2.len());
        for (p1, p2) in points1.iter().zip(points2.iter()) {
            assert!((p1.x - p2.x).abs() < 1e-10);
            assert!((p1.y - p2.y).abs() < 1e-10);
        }
    }

    #[test]
    fn test_lloyd_relaxation() {
        let config = PointGeneratorConfig {
            num_points: 50,
            lloyd_iterations: 2,
            ..Default::default()
        };
        let points = generate_relaxed_points(&config);
        assert_eq!(points.len(), 50);

        // All points should still be within bounds
        for p in &points {
            assert!(config.bounds.contains(p));
        }
    }
}
