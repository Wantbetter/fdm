extern crate fdm;

use ndarray::{Array2, Array, ShapeBuilder};
use fdm::sw;

fn main() {
    let v = sw::cs(5);
    dbg!(&v);

    let r = sw::is_stable(1000.0, 0.00025, 1.0, 1.0, &v);
    dbg!(r);
}
