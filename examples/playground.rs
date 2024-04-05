use ddtri::formula::{distance_difference_triangulation, BeaconDistances, DistanceDifferences};
use glam::{dvec2, DVec2};

fn main() {
    let beacon0 = dvec2(0., 0.);
    let beacon1 = dvec2(7., 0.);
    let beacon2 = dvec2(5., 6.);

    let position = dvec2(4., 2.) + (2. * rand::random::<DVec2>() - dvec2(1., 1.));

    let dist0 = beacon0.distance(position);
    let dist1 = beacon1.distance(position);
    let dist2 = beacon2.distance(position);

    let distdiffs = DistanceDifferences {
        ddto0to1: dist0 - dist1,
        ddto0to2: dist0 - dist2,
    };

    let beacondists = BeaconDistances {
        d0to1: beacon0.distance(beacon1),
        d0to2: beacon0.distance(beacon2),
        d1to2: beacon1.distance(beacon2),
    };

    let res: DVec2 = distance_difference_triangulation(distdiffs, beacondists).into();

    println!("True Pos:    {}", position);
    println!("Inferred Pos:{}", res);
}
