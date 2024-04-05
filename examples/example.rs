use ddtri::{distance_difference_triangulation, BeaconDistances, DistanceDifferences};
use rand::{
    distributions::{Distribution, Uniform},
    thread_rng,
};

fn main() {
    let beacon0 = (0., 0.);
    let beacon1 = (7., 0.);
    let beacon2 = (5., 6.);

    let position = (
        4. + Uniform::from(-1.0..1.0).sample(&mut thread_rng()),
        2. + Uniform::from(-1.0..1.0).sample(&mut thread_rng()),
    );

    let dist0 = distance(beacon0, position);
    let dist1 = distance(beacon1, position);
    let dist2 = distance(beacon2, position);

    let position_estimate = distance_difference_triangulation(
        DistanceDifferences {
            ddto0to1: dist0 - dist1,
            ddto0to2: dist0 - dist2,
        },
        BeaconDistances {
            d0to1: distance(beacon0, beacon1),
            d0to2: distance(beacon0, beacon2),
            d1to2: distance(beacon1, beacon2),
        },
    );

    println!("Original Position:");
    println!("  x: {}", position.0);
    println!("  y: {}", position.1);

    println!("Inferred Position:");
    println!("  x: {}", position_estimate.0);
    println!("  y: {}", position_estimate.1);

    let err = (
        (position_estimate.0 - position.0).abs(),
        (position_estimate.1 - position.1).abs(),
    );
    let err_total = distance((0., 0.), err);

    println!("Error:");
    println!("  x:     {}", err.0);
    println!("  y:     {}", err.1);
    println!("  total: {}", err_total);
}

/// Euclidian distance between two points.
fn distance(p: (f64, f64), q: (f64, f64)) -> f64 {
    let d0 = p.0 - q.0;
    let d1 = p.1 - q.1;
    (d0 * d0 + d1 * d1).sqrt()
}
