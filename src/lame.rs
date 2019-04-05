use ndarray::{Array, Array2, ShapeBuilder};

// #[derive(Debug, Clone)]
// pub struct Lame(Array2<f64>, Array2<f64>);

// Lame.0 => lambda
// Lame.1 => miu
pub fn lame(vp: &Array2<f64>, vs: &Array2<f64>, p: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
    assert!(vp.shape() == vs.shape() && vs.shape() == p.shape());

    let rows = vp.rows();
    let cols = vp.cols();

    let lame1 = Array2::<f64>::from_shape_fn((rows, cols).f(), |ix| {
        let j = ix.0;
        let i = ix.1;
        p[[j, i]] * (vp[[j, i]].powf(2.0) - 2.0 * vs[[j, i]].powf(2.0))
    });

    let lame2 = Array2::<f64>::from_shape_fn((rows, cols).f(), |ix| {
        let j = ix.0;
        let i = ix.1;
        p[[j, i]] * vs[[j, i]].powf(2.0)
    });

    (lame1, lame2)
}

// #[test]
// pub fn test_lame() {
//     let p = Array::from_elem((5, 6).f(), 2000.0);
//     let vp = Array::from_elem((5, 6).f(), 1000.0);
//     let vs = Array::from_elem((5, 6).f(), 570.0);

//     let lame = lame(&vp, &vs, &p);

//     assert_eq!(lame.0[[0, 0]], 700400000.0f64);
//     assert_eq!(lame.1[[0, 0]], 649800000.0f64);
// }
