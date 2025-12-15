//! Basic geometric types and utilities.

use std::ops::{Add, Sub, Mul, Div};

/// A 2D point/vector.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub const ZERO: Point = Point { x: 0.0, y: 0.0 };

    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Squared distance to another point.
    #[inline]
    pub fn distance_sq(&self, other: &Point) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    /// Distance to another point.
    #[inline]
    pub fn distance(&self, other: &Point) -> f64 {
        self.distance_sq(other).sqrt()
    }

    /// Length of the vector.
    #[inline]
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Normalized vector.
    #[inline]
    pub fn normalize(&self) -> Self {
        let len = self.length();
        if len > 1e-10 {
            Self { x: self.x / len, y: self.y / len }
        } else {
            Self::ZERO
        }
    }

    /// Linear interpolation.
    #[inline]
    pub fn lerp(&self, other: &Point, t: f64) -> Self {
        Self {
            x: self.x + (other.x - self.x) * t,
            y: self.y + (other.y - self.y) * t,
        }
    }
}

impl Add for Point {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self { x: self.x + rhs.x, y: self.y + rhs.y }
    }
}

impl Sub for Point {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self { x: self.x - rhs.x, y: self.y - rhs.y }
    }
}

impl Mul<f64> for Point {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        Self { x: self.x * rhs, y: self.y * rhs }
    }
}

impl Div<f64> for Point {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f64) -> Self {
        Self { x: self.x / rhs, y: self.y / rhs }
    }
}

impl From<delaunator::Point> for Point {
    fn from(p: delaunator::Point) -> Self {
        Self { x: p.x, y: p.y }
    }
}

impl From<Point> for delaunator::Point {
    fn from(p: Point) -> Self {
        Self { x: p.x, y: p.y }
    }
}

/// Compute the circumcenter of a triangle.
/// The circumcenter is equidistant from all three vertices.
pub fn circumcenter(a: &Point, b: &Point, c: &Point) -> Point {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let ex = c.x - a.x;
    let ey = c.y - a.y;

    let bl = dx * dx + dy * dy;
    let cl = ex * ex + ey * ey;
    let d = dx * ey - dy * ex;

    if d.abs() < 1e-10 {
        // Degenerate triangle, return centroid instead
        return Point::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
        );
    }

    let d = 0.5 / d;
    let x = a.x + (ey * bl - dy * cl) * d;
    let y = a.y + (dx * cl - ex * bl) * d;

    Point::new(x, y)
}

/// Compute the centroid of a polygon given its vertices.
pub fn polygon_centroid(vertices: &[Point]) -> Point {
    if vertices.is_empty() {
        return Point::ZERO;
    }
    if vertices.len() == 1 {
        return vertices[0];
    }
    if vertices.len() == 2 {
        return vertices[0].lerp(&vertices[1], 0.5);
    }

    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut area = 0.0;

    for i in 0..vertices.len() {
        let j = (i + 1) % vertices.len();
        let cross = vertices[i].x * vertices[j].y - vertices[j].x * vertices[i].y;
        area += cross;
        cx += (vertices[i].x + vertices[j].x) * cross;
        cy += (vertices[i].y + vertices[j].y) * cross;
    }

    area *= 0.5;
    if area.abs() < 1e-10 {
        // Degenerate polygon, return simple average
        let sum: Point = vertices.iter().fold(Point::ZERO, |acc, p| acc + *p);
        return sum / vertices.len() as f64;
    }

    let factor = 1.0 / (6.0 * area);
    Point::new(cx * factor, cy * factor)
}

/// Axis-aligned bounding box.
#[derive(Debug, Clone, Copy)]
pub struct BoundingBox {
    pub min: Point,
    pub max: Point,
}

impl BoundingBox {
    pub fn new(min_x: f64, min_y: f64, max_x: f64, max_y: f64) -> Self {
        Self {
            min: Point::new(min_x, min_y),
            max: Point::new(max_x, max_y),
        }
    }

    pub fn width(&self) -> f64 {
        self.max.x - self.min.x
    }

    pub fn height(&self) -> f64 {
        self.max.y - self.min.y
    }

    pub fn center(&self) -> Point {
        Point::new(
            (self.min.x + self.max.x) / 2.0,
            (self.min.y + self.max.y) / 2.0,
        )
    }

    pub fn contains(&self, p: &Point) -> bool {
        p.x >= self.min.x && p.x <= self.max.x &&
        p.y >= self.min.y && p.y <= self.max.y
    }

    /// Check if point is on the boundary (within epsilon).
    pub fn is_on_boundary(&self, p: &Point, epsilon: f64) -> bool {
        let on_left = (p.x - self.min.x).abs() < epsilon;
        let on_right = (p.x - self.max.x).abs() < epsilon;
        let on_bottom = (p.y - self.min.y).abs() < epsilon;
        let on_top = (p.y - self.max.y).abs() < epsilon;
        on_left || on_right || on_bottom || on_top
    }
}
