pub struct DistanceDifferences {
    pub ddto0to1: f64,
    pub ddto0to2: f64,
}
pub struct BeaconDistances {
    pub d0to1: f64,
    pub d0to2: f64,
    pub d1to2: f64,
}

pub fn distance_difference_triangulation(
    distance_differences: DistanceDifferences,
    beacon_distances: BeaconDistances,
) -> (f64, f64) {
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
    let x2 = (bd02 * bd02 - bd12 * bd12 + bd01 * bd01) / (2. * bd01);
    let y2 = (bd02 * bd02 - x2 * x2).sqrt();

    bigboy(x1, x2, y2, a, b)
}

/// See https://github.com/Garbaz/triangulation_from_dist_diff/. Note that the
/// indices are shifted down by one here for consistency of notation.
fn bigboy(x1: f64, x2: f64, y2: f64, a: f64, b: f64) -> (f64, f64) {
    let q = 4. * x1 * x1 * a * a + 4. * x2 * x2 * a * a + 4. * y2 * y2 * a * a
        - 8. * x1 * x2 * a * a
        + 8. * b * x1 * x1 * a
        - 8. * b * x1 * x2 * a
        + 4. * b * b * x1 * x1
        - 4. * x1 * x1 * y2 * y2;

    let x = -((4. * x1 * x1 * a * a * a * a) / q
        + (4. * x2 * x2 * a * a * a * a) / q
        + (4. * y2 * y2 * a * a * a * a) / q
        - (8. * x1 * x2 * a * a * a * a) / q
        + (12. * b * x1 * x1 * a * a * a) / q
        + (8. * b * x2 * x2 * a * a * a) / q
        + (8. * b * y2 * y2 * a * a * a) / q
        - (20. * b * x1 * x2 * a * a * a) / q
        - (4. * x1 * x2 * x2 * x2 * a * a) / q
        + (12. * b * b * x1 * x1 * a * a) / q
        + (8. * x1 * x1 * x2 * x2 * a * a) / q
        - (4. * x1 * x2 * y2 * y2 * a * a) / q
        - (4. * x1 * x1 * x1 * x2 * a * a) / q
        - (12. * b * b * x1 * x2 * a * a) / q
        - a * a
        - 2. * b * a
        + (4. * b * b * b * x1 * x1 * a) / q
        + (4. * b * x1 * x1 * x2 * x2 * a) / q
        - (4. * b * x1 * x1 * y2 * y2 * a) / q
        - (4. * b * x1 * x1 * x1 * x2 * a) / q
        - x1 * x1
        + (4.
            * (y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                - a * a * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                - b * b * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                + x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                - 2. * a * b * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
                - 2. * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                - 2. * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                + 2. * a * a * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                + 2. * b * b * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                + 4. * a * b * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
                + y2 * y2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                - 2. * a * a * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                - 2. * b * b * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                + 2. * x2 * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                - 2. * a * b * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
                + a * a * a * a * y2 * y2 * x1 * x1 * x1 * x1
                + b * b * b * b * y2 * y2 * x1 * x1 * x1 * x1
                + x2 * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
                + 2. * a * b * b * b * y2 * y2 * x1 * x1 * x1 * x1
                + 2. * a * a * b * b * y2 * y2 * x1 * x1 * x1 * x1
                - 2. * a * a * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
                - 2. * b * b * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
                - 2. * a * b * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
                + 2. * a * a * a * b * y2 * y2 * x1 * x1 * x1 * x1
                + 2. * a * a * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1
                + 2. * a * a * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1
                - 2. * a * a * a * a * x2 * y2 * y2 * x1 * x1 * x1
                - 2. * a * a * b * b * x2 * y2 * y2 * x1 * x1 * x1
                - 4. * a * a * a * b * x2 * y2 * y2 * x1 * x1 * x1
                - a * a * y2 * y2 * y2 * y2 * y2 * y2 * x1 * x1
                + a * a * a * a * y2 * y2 * y2 * y2 * x1 * x1
                + 2. * a * a * b * b * y2 * y2 * y2 * y2 * x1 * x1
                - 2. * a * a * x2 * x2 * y2 * y2 * y2 * y2 * x1 * x1
                + 2. * a * a * a * b * y2 * y2 * y2 * y2 * x1 * x1
                - a * a * b * b * b * b * y2 * y2 * x1 * x1
                - a * a * x2 * x2 * x2 * x2 * y2 * y2 * x1 * x1
                - 2. * a * a * a * b * b * b * y2 * y2 * x1 * x1
                - a * a * a * a * b * b * y2 * y2 * x1 * x1
                + a * a * a * a * x2 * x2 * y2 * y2 * x1 * x1
                + 2. * a * a * b * b * x2 * x2 * y2 * y2 * x1 * x1
                + 2. * a * a * a * b * x2 * x2 * y2 * y2 * x1 * x1)
                .sqrt()
            * a)
            / q)
        / (2. * x1);

    let qq = (y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        - a * a * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        - b * b * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        + x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        - 2. * a * b * y2 * y2 * x1 * x1 * x1 * x1 * x1 * x1
        - 2. * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        - 2. * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        + 2. * a * a * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        + 2. * b * b * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        + 4. * a * b * x2 * y2 * y2 * x1 * x1 * x1 * x1 * x1
        + y2 * y2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        - 2. * a * a * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        - 2. * b * b * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        + 2. * x2 * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        - 2. * a * b * y2 * y2 * y2 * y2 * x1 * x1 * x1 * x1
        + a * a * a * a * y2 * y2 * x1 * x1 * x1 * x1
        + b * b * b * b * y2 * y2 * x1 * x1 * x1 * x1
        + x2 * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
        + 2. * a * b * b * b * y2 * y2 * x1 * x1 * x1 * x1
        + 2. * a * a * b * b * y2 * y2 * x1 * x1 * x1 * x1
        - 2. * a * a * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
        - 2. * b * b * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
        - 2. * a * b * x2 * x2 * y2 * y2 * x1 * x1 * x1 * x1
        + 2. * a * a * a * b * y2 * y2 * x1 * x1 * x1 * x1
        + 2. * a * a * x2 * y2 * y2 * y2 * y2 * x1 * x1 * x1
        + 2. * a * a * x2 * x2 * x2 * y2 * y2 * x1 * x1 * x1
        - 2. * a * a * a * a * x2 * y2 * y2 * x1 * x1 * x1
        - 2. * a * a * b * b * x2 * y2 * y2 * x1 * x1 * x1
        - 4. * a * a * a * b * x2 * y2 * y2 * x1 * x1 * x1
        - a * a * y2 * y2 * y2 * y2 * y2 * y2 * x1 * x1
        + a * a * a * a * y2 * y2 * y2 * y2 * x1 * x1
        + 2. * a * a * b * b * y2 * y2 * y2 * y2 * x1 * x1
        - 2. * a * a * x2 * x2 * y2 * y2 * y2 * y2 * x1 * x1
        + 2. * a * a * a * b * y2 * y2 * y2 * y2 * x1 * x1
        - a * a * b * b * b * b * y2 * y2 * x1 * x1
        - a * a * x2 * x2 * x2 * x2 * y2 * y2 * x1 * x1
        - 2. * a * a * a * b * b * b * y2 * y2 * x1 * x1
        - a * a * a * a * b * b * y2 * y2 * x1 * x1
        + a * a * a * a * x2 * x2 * y2 * y2 * x1 * x1
        + 2. * a * a * b * b * x2 * x2 * y2 * y2 * x1 * x1
        + 2. * a * a * a * b * x2 * x2 * y2 * y2 * x1 * x1)
        .sqrt();

    let y = -(-x1 * a * a + x2 * a * a - 2. * b * x1 * a + 2. * b * x2 * a
        - (x1
            * (-4. * x1 * x1 * a * a * a - 4. * x2 * x2 * a * a * a - 4. * y2 * y2 * a * a * a
                + 8. * x1 * x2 * a * a * a
                - 12. * b * x1 * x1 * a * a
                - 8. * b * x2 * x2 * a * a
                - 8. * b * y2 * y2 * a * a
                + 20. * b * x1 * x2 * a * a
                + 4. * x1 * x2 * x2 * x2 * a
                - 12. * b * b * x1 * x1 * a
                - 8. * x1 * x1 * x2 * x2 * a
                + 4. * x1 * x2 * y2 * y2 * a
                + 4. * x1 * x1 * x1 * x2 * a
                + 12. * b * b * x1 * x2 * a
                - 4. * b * b * b * x1 * x1
                - 4. * b * x1 * x1 * x2 * x2
                + 4. * b * x1 * x1 * y2 * y2
                + 4. * b * x1 * x1 * x1 * x2
                - 4. * qq)
            * a)
            / q
        + (x2
            * (-4. * x1 * x1 * a * a * a - 4. * x2 * x2 * a * a * a - 4. * y2 * y2 * a * a * a
                + 8. * x1 * x2 * a * a * a
                - 12. * b * x1 * x1 * a * a
                - 8. * b * x2 * x2 * a * a
                - 8. * b * y2 * y2 * a * a
                + 20. * b * x1 * x2 * a * a
                + 4. * x1 * x2 * x2 * x2 * a
                - 12. * b * b * x1 * x1 * a
                - 8. * x1 * x1 * x2 * x2 * a
                + 4. * x1 * x2 * y2 * y2 * a
                + 4. * x1 * x1 * x1 * x2 * a
                + 12. * b * b * x1 * x2 * a
                - 4. * b * b * b * x1 * x1
                - 4. * b * x1 * x1 * x2 * x2
                + 4. * b * x1 * x1 * y2 * y2
                + 4. * b * x1 * x1 * x1 * x2
                - 4. * qq)
            * a)
            / q
        - x1 * x2 * x2
        - x1 * y2 * y2
        - b * b * x1
        + x1 * x1 * x2
        - (b * x1
            * (-4. * x1 * x1 * a * a * a - 4. * x2 * x2 * a * a * a - 4. * y2 * y2 * a * a * a
                + 8. * x1 * x2 * a * a * a
                - 12. * b * x1 * x1 * a * a
                - 8. * b * x2 * x2 * a * a
                - 8. * b * y2 * y2 * a * a
                + 20. * b * x1 * x2 * a * a
                + 4. * x1 * x2 * x2 * x2 * a
                - 12. * b * b * x1 * x1 * a
                - 8. * x1 * x1 * x2 * x2 * a
                + 4. * x1 * x2 * y2 * y2 * a
                + 4. * x1 * x1 * x1 * x2 * a
                + 12. * b * b * x1 * x2 * a
                - 4. * b * b * b * x1 * x1
                - 4. * b * x1 * x1 * x2 * x2
                + 4. * b * x1 * x1 * y2 * y2
                + 4. * b * x1 * x1 * x1 * x2
                - 4. * qq))
            / q)
        / (2. * x1 * y2);

    (x, y)
}
