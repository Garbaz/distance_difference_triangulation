# Distance Difference Triangulation

‖ [__Docs.rs__](https://docs.rs/ddtri) ‖ [__Lib.rs__](https://lib.rs/crates/ddtri) ‖ [__Crates.io__](https://crates.io/crates/ddtri/) ‖


This crate consists of exactly one function: [`distance_difference_triangulation`].

The premise is as follows: You are at some unknown position in 2D space. There
are three beacons. You _do not_ know the distances to these three beacons,
otherwise you could just do normal triangulation and be done with it. But you
_do_ know the _differences between the distances_ to these beacons. You also
know the distances between the beacons themselves. From this information, the
function this crate provides computes your position, relative to the coordinate
system defined by the beacons.

In mathematical terms:

Your _unknown_ position is `(x,y)`.

The _unknown_ distances to the three beacons are `d0`, `d1` and `d2`.

You do know `dd01 = d0 - d1` and `dd02 = d0 - d2` (and also `dd12 = d1 - d2`, but that's redundant).

You also know the distances between the beacons `d01`, `d02` and `d12`.

The function of this crate computes for you your position _in_ the coordinate
system where beacon 0 is at `(0,0)` and beacon 1 is on the x axis.

See also [triangulation_from_dist_diff](https://github.com/Garbaz/triangulation_from_dist_diff/).
