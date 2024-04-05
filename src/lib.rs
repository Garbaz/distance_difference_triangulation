#![doc = include_str!("../README.md")]

use std::ops::{Add, Div, Mul, Neg, Sub};

/// The distance differences.
///
/// If the point we want to locate is `(x,y)`, and the distance from our point
/// `(x,y)` to beacon0 is `d0`, to beacon1 `d1`, and to beacon2 `d2`, then:
///
/// `ddto0to1 == d0 - d1` and `ddto0to2 == d0 - d2`.
pub struct DistanceDifferences<F> {
    pub ddto0to1: F,
    pub ddto0to2: F,
}

/// The distances between the beacons.
///
/// `d0to1` is the distance between beacon 0 and beacon 1.\
/// `d0to2` is the distance between beacon 0 and beacon 2.\
/// `d1to2` is the distance between beacon 1 and beacon 2.
pub struct BeaconDistances<F> {
    pub d0to1: F,
    pub d0to2: F,
    pub d1to2: F,
}

/// Triangulate the position of an unknown point `(x,y)` from the distance
/// differences to three beacons and the distances between the three beacons.
pub fn distance_difference_triangulation<F>(
    distance_differences: DistanceDifferences<F>,
    beacon_distances: BeaconDistances<F>,
) -> (F, F)
where
    F: Sqrt
        + Neg<Output = F>
        + Add<Output = F>
        + Sub<Output = F>
        + Mul<Output = F>
        + Div<Output = F>
        + From<u8>
        + Copy,
{
    let DistanceDifferences {
        ddto0to1: dd01,
        ddto0to2: dd02,
    } = distance_differences;
    let BeaconDistances {
        d0to1: bd01,
        d0to2: bd02,
        d1to2: bd12,
    } = beacon_distances;

    let a = dd01;
    let b = dd02 - dd01;

    let x1 = bd01;
    let x2 = (bd02 * bd02 - bd12 * bd12 + bd01 * bd01) / (F::from(2) * bd01);
    let y2 = (bd02 * bd02 - x2 * x2)._sqrt();

    bigboy(x1, x2, y2, a, b)
}

/// See https://github.com/Garbaz/triangulation_from_dist_diff/. Note that the
/// indices are shifted down by one here for consistency of notation.
fn bigboy<F>(x1: F, x2: F, y2: F, a: F, b: F) -> (F, F)
where
    F: Sqrt
        + Neg<Output = F>
        + Add<Output = F>
        + Sub<Output = F>
        + Mul<Output = F>
        + Div<Output = F>
        + From<u8>
        + Copy,
{
    let q =
        F::from(4) * x1 * x1 * a * a + F::from(4) * x2 * x2 * a * a + F::from(4) * y2 * y2 * a * a
            - F::from(8) * x1 * x2 * a * a
            + F::from(8) * b * x1 * x1 * a
            - F::from(8) * b * x1 * x2 * a
            + F::from(4) * b * b * x1 * x1
            - F::from(4) * x1 * x1 * y2 * y2;

    let x = -((F::from(4) * x1 * x1 * a * a * a * a) / q
        + (F::from(4) * x2 * x2 * a * a * a * a) / q
        + (F::from(4) * y2 * y2 * a * a * a * a) / q
        - (F::from(8) * x1 * x2 * a * a * a * a) / q
        + (F::from(12) * b * x1 * x1 * a * a * a) / q
        + (F::from(8) * b * x2 * x2 * a * a * a) / q
        + (F::from(8) * b * y2 * y2 * a * a * a) / q
        - (F::from(20) * b * x1 * x2 * a * a * a) / q
        - (F::from(4) * x1 * x2 * x2 * x2 * a * a) / q
        + (F::from(12) * b * b * x1 * x1 * a * a) / q
        + (F::from(8) * x1 * x1 * x2 * x2 * a * a) / q
        - (F::from(4) * x1 * x2 * y2 * y2 * a * a) / q
        - (F::from(4) * x1 * x1 * x1 * x2 * a * a) / q
        - (F::from(12) * b * b * x1 * x2 * a * a) / q
        - a * a
        - F::from(2) * b * a
        + (F::from(4) * b * b * b * x1 * x1 * a) / q
        + (F::from(4) * b * x1 * x1 * x2 * x2 * a) / q
        - (F::from(4) * b * x1 * x1 * y2 * y2 * a) / q
        - (F::from(4) * b * x1 * x1 * x1 * x2 * a) / q
        - x1 * x1
        + (F::from(4)
            * (y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                - a * a * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                - b * b * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                + x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                - F::from(2) * a * b * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                - F::from(2) * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                - F::from(2) * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                + F::from(2) * a * a * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                + F::from(2) * b * b * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                + F::from(4) * a * b * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                + y2 * y2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                - F::from(2) * a * a * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                - F::from(2) * b * b * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                + F::from(2) * x2 * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                - F::from(2) * a * b * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                + a * a * a * a * y2 * y2 * x1 * x1 * x1 * x1
                + b * b * b * b * y2 * y2 * x1 * x1 * x1 * x1
                + x2 * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
                + F::from(2) * a * b * b * b * y2 * y2 * x1 * x1 * x1 * x1
                + F::from(2) * a * a * b * b * y2 * y2 * x1 * x1 * x1 * x1
                - F::from(2) * a * a * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
                - F::from(2) * b * b * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
                - F::from(2) * a * b * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
                + F::from(2) * a * a * a * b * y2 * y2 * x1 * x1 * x1 * x1
                + F::from(2) * a * a * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1
                + F::from(2) * a * a * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1
                - F::from(2) * a * a * a * a * x2 * y2 * y2 * x1 * x1 * x1
                - F::from(2) * a * a * b * b * x2 * y2 * y2 * x1 * x1 * x1
                - F::from(4) * a * a * a * b * x2 * y2 * y2 * x1 * x1 * x1
                - a * a * y2 * y2 * y2 * y2 * y2 * y2 * x1 * x1
                + a * a * a * a * y2 * y2 * y2 * y2 * x1 * x1
                + F::from(2) * a * a * b * b * y2 * y2 * y2 * y2 * x1 * x1
                - F::from(2) * a * a * x2 * x2 * y2 * y2 * y2 * y2 * x1 * x1
                + F::from(2) * a * a * a * b * y2 * y2 * y2 * y2 * x1 * x1
                - a * a * b * b * b * b * y2 * y2 * x1 * x1
                - a * a * x2 * x2 * x2 * x2 * y2 * y2 * x1 * x1
                - F::from(2) * a * a * a * b * b * b * y2 * y2 * x1 * x1
                - a * a * a * a * b * b * y2 * y2 * x1 * x1
                + a * a * a * a * x2 * x2 * y2 * y2 * x1 * x1
                + F::from(2) * a * a * b * b * x2 * x2 * y2 * y2 * x1 * x1
                + F::from(2) * a * a * a * b * x2 * x2 * y2 * y2 * x1 * x1)
                ._sqrt()
            * a)
            / q)
        / (F::from(2) * x1);

    let qq = (y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        - a * a * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        - b * b * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        + x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        - F::from(2) * a * b * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        - F::from(2) * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        - F::from(2) * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        + F::from(2) * a * a * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        + F::from(2) * b * b * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        + F::from(4) * a * b * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        + y2 * y2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        - F::from(2) * a * a * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        - F::from(2) * b * b * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        + F::from(2) * x2 * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        - F::from(2) * a * b * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        + a * a * a * a * y2 * y2 * x1 * x1 * x1 * x1
        + b * b * b * b * y2 * y2 * x1 * x1 * x1 * x1
        + x2 * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
        + F::from(2) * a * b * b * b * y2 * y2 * x1 * x1 * x1 * x1
        + F::from(2) * a * a * b * b * y2 * y2 * x1 * x1 * x1 * x1
        - F::from(2) * a * a * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
        - F::from(2) * b * b * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
        - F::from(2) * a * b * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
        + F::from(2) * a * a * a * b * y2 * y2 * x1 * x1 * x1 * x1
        + F::from(2) * a * a * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1
        + F::from(2) * a * a * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1
        - F::from(2) * a * a * a * a * x2 * y2 * y2 * x1 * x1 * x1
        - F::from(2) * a * a * b * b * x2 * y2 * y2 * x1 * x1 * x1
        - F::from(4) * a * a * a * b * x2 * y2 * y2 * x1 * x1 * x1
        - a * a * y2 * y2 * y2 * y2 * y2 * y2 * x1 * x1
        + a * a * a * a * y2 * y2 * y2 * y2 * x1 * x1
        + F::from(2) * a * a * b * b * y2 * y2 * y2 * y2 * x1 * x1
        - F::from(2) * a * a * x2 * x2 * y2 * y2 * y2 * y2 * x1 * x1
        + F::from(2) * a * a * a * b * y2 * y2 * y2 * y2 * x1 * x1
        - a * a * b * b * b * b * y2 * y2 * x1 * x1
        - a * a * x2 * x2 * x2 * x2 * y2 * y2 * x1 * x1
        - F::from(2) * a * a * a * b * b * b * y2 * y2 * x1 * x1
        - a * a * a * a * b * b * y2 * y2 * x1 * x1
        + a * a * a * a * x2 * x2 * y2 * y2 * x1 * x1
        + F::from(2) * a * a * b * b * x2 * x2 * y2 * y2 * x1 * x1
        + F::from(2) * a * a * a * b * x2 * x2 * y2 * y2 * x1 * x1)
        ._sqrt();

    let y = -(-x1 * a * a + x2 * a * a - F::from(2) * b * x1 * a + F::from(2) * b * x2 * a
        - x1 * (-(F::from(4)) * x1 * x1 * a * a * a
            - F::from(4) * x2 * x2 * a * a * a
            - F::from(4) * y2 * y2 * a * a * a
            + F::from(8) * x1 * x2 * a * a * a
            - F::from(12) * b * x1 * x1 * a * a
            - F::from(8) * b * x2 * x2 * a * a
            - F::from(8) * b * y2 * y2 * a * a
            + F::from(20) * b * x1 * x2 * a * a
            + F::from(4) * x1 * x2 * x2 * x2 * a
            - F::from(12) * b * b * x1 * x1 * a
            - F::from(8) * x1 * x1 * x2 * x2 * a
            + F::from(4) * x1 * x2 * y2 * y2 * a
            + F::from(4) * x1 * x1 * x1 * x2 * a
            + F::from(12) * b * b * x1 * x2 * a
            - F::from(4) * b * b * b * x1 * x1
            - F::from(4) * b * x1 * x1 * x2 * x2
            + F::from(4) * b * x1 * x1 * y2 * y2
            + F::from(4) * b * x1 * x1 * x1 * x2
            - F::from(4) * qq)
            * a
            / q
        + (x2
            * (-(F::from(4)) * x1 * x1 * a * a * a
                - F::from(4) * x2 * x2 * a * a * a
                - F::from(4) * y2 * y2 * a * a * a
                + F::from(8) * x1 * x2 * a * a * a
                - F::from(12) * b * x1 * x1 * a * a
                - F::from(8) * b * x2 * x2 * a * a
                - F::from(8) * b * y2 * y2 * a * a
                + F::from(20) * b * x1 * x2 * a * a
                + F::from(4) * x1 * x2 * x2 * x2 * a
                - F::from(12) * b * b * x1 * x1 * a
                - F::from(8) * x1 * x1 * x2 * x2 * a
                + F::from(4) * x1 * x2 * y2 * y2 * a
                + F::from(4) * x1 * x1 * x1 * x2 * a
                + F::from(12) * b * b * x1 * x2 * a
                - F::from(4) * b * b * b * x1 * x1
                - F::from(4) * b * x1 * x1 * x2 * x2
                + F::from(4) * b * x1 * x1 * y2 * y2
                + F::from(4) * b * x1 * x1 * x1 * x2
                - F::from(4) * qq)
            * a)
            / q
        - x1 * x2 * x2
        - x1 * y2 * y2
        - b * b * x1
        + x1 * x1 * x2
        - (b * x1
            * (-(F::from(4)) * x1 * x1 * a * a * a
                - F::from(4) * x2 * x2 * a * a * a
                - F::from(4) * y2 * y2 * a * a * a
                + F::from(8) * x1 * x2 * a * a * a
                - F::from(12) * b * x1 * x1 * a * a
                - F::from(8) * b * x2 * x2 * a * a
                - F::from(8) * b * y2 * y2 * a * a
                + F::from(20) * b * x1 * x2 * a * a
                + F::from(4) * x1 * x2 * x2 * x2 * a
                - F::from(12) * b * b * x1 * x1 * a
                - F::from(8) * x1 * x1 * x2 * x2 * a
                + F::from(4) * x1 * x2 * y2 * y2 * a
                + F::from(4) * x1 * x1 * x1 * x2 * a
                + F::from(12) * b * b * x1 * x2 * a
                - F::from(4) * b * b * b * x1 * x1
                - F::from(4) * b * x1 * x1 * x2 * x2
                + F::from(4) * b * x1 * x1 * y2 * y2
                + F::from(4) * b * x1 * x1 * x1 * x2
                - F::from(4) * qq))
            / q)
        / (F::from(2) * x1 * y2);

    (x, y)
}

/// A trait for types that implement the the square-root function.
pub trait Sqrt {
    /// Returns the square root of a number.
    ///
    /// Returns NaN if `self` is a negative number other than `-0.0`.
    fn _sqrt(self) -> Self;
}

impl Sqrt for f32 {
    fn _sqrt(self) -> Self {
        f32::sqrt(self)
    }
}

impl Sqrt for f64 {
    fn _sqrt(self) -> Self {
        f64::sqrt(self)
    }
}
