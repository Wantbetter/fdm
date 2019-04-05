extern crate fdm;
#[macro_use(s)]
extern crate ndarray;

use ndarray::{Array2, Array, ShapeBuilder, Axis};
use fdm::sw;

use std::ops::{Div, Mul};

fn main() {
    let r = 1e-6f64;

    let vsmax = 1000.0;
    let nx = 159;
    let nz = 159;
    let pml_h = 20;
    let nz_pml = nz + 2 * pml_h;
    let nx_pml = nx + 2 * pml_h;
    let dx = 5.0;
    let dz = 5.0;
    let dt = 1e-3;
    let vp = Array2::<f64>::from_shape_fn((nz_pml, nx_pml).f(), |ix| 2000.0);
    let vs = Array2::<f64>::from_shape_fn((nz_pml, nx_pml).f(), |ix| 1000.0);
    let pp = Array2::<f64>::from_shape_fn((nz_pml, nx_pml).f(), |ix| 1000.0);

    let d_pml_x0 = (1.0 / r).ln() * 3.0 * vsmax / (2.0 * pml_h as f64);
    
    let mut d_pml_x = Array2::<f64>::zeros((nz_pml, nx_pml));

    let mat1 = Array2::<f64>::ones((nz_pml, 1).f());
    let mat2 = Array2::<f64>::from_shape_fn((1, pml_h).f(), |ix| {
        let j = ix.1;
        ((pml_h - j) as f64 / pml_h as f64).powf(2.0)
    });

    d_pml_x.slice_mut(s![.., 0..pml_h]).assign(&mat1.dot(&(d_pml_x0 * mat2)));
    let d_pml_x_tmp = d_pml_x.slice(s![.., 0..pml_h; -1]).to_owned();  //解决
    d_pml_x.slice_mut(s![.., nx + pml_h..nx_pml]).assign(&d_pml_x_tmp); 

    let d_pml_z0 = d_pml_x0;
    let mut d_pml_z = Array2::<f64>::zeros((nz_pml, nx_pml));
    let mat1 = Array2::<f64>::ones((1, nx_pml).f());
    let mat2 = Array2::<f64>::from_shape_fn((pml_h, 1).f(), |ix| {
        let j = ix.0;
        ((pml_h - j) as f64 / pml_h as f64).powf(2.0)
    });

    d_pml_z.slice_mut(s![0..pml_h, ..]).assign(&((d_pml_z0 * mat2).dot(&mat1)));
    let d_pml_z_tmp = d_pml_z.slice(s![0..pml_h;-1, ..]).to_owned();  //解决
    d_pml_z.slice_mut(s![nz + pml_h..nz_pml, ..]).assign(&d_pml_z_tmp);

    let (lambda, miu) = fdm::lame::lame(&vp, &vs, &pp);

    let lame_c = lambda.clone() + 2.0 * miu.clone();

    let c_t = (2.0 * dt).div(2.0 + dt * &d_pml_x);

    dbg!(lame_c.shape());
    dbg!(c_t.shape());
    
    let c1 = (2.0 - dt * &d_pml_x).div(2.0 + dt * &d_pml_x);
    let c2 = (2.0 - dt * &d_pml_z).div(2.0 + dt * &d_pml_z);
    let c3 = (2.0 * dt).div(2.0 + dt * &d_pml_x).div(&pp).div(dx);
    let c4 = (2.0 * dt).div(2.0 + dt * &d_pml_z).div(&pp).div(dz);
    let c5 = (2.0 * dt).div(2.0 + dt * &d_pml_x) * (&lame_c).div(dx);
    let c6 = (2.0 * dt).div(2.0 + dt * &d_pml_z) * (&lambda).div(dz);
    let c7 = (2.0 * dt).div(2.0 + dt * &d_pml_x) * (&lambda).div(dx);
    let c8 = (2.0 * dt).div(2.0 + dt * &d_pml_z) * (&lame_c).div(dz);
    let c9 = (2.0 * dt).div(2.0 + dt * &d_pml_x) * (&miu).div(dx);
    let c0 = (2.0 * dt).div(2.0 + dt * d_pml_z) * (&miu).div(dz);
    // dbg!(4000000000.0f64 * c_t.slice(s![0, ..]).to_owned().sum_axis(Axis(0)));
    dbg!(c5.slice(s![0, ..]));
    dbg!(c6.slice(s![0, ..]));
    dbg!(c7.slice(s![0, ..]));
    dbg!(c8.slice(s![0, ..]));
    dbg!(c9.slice(s![0, ..]));
    dbg!(c0.slice(s![0, ..]));
}
