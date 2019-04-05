use ndarray::{Array, Array2, Array3, ShapeBuilder};
use std::{ops::Div, path::Path};

use crate::io;

pub fn cs(jl: usize) -> Vec<f64> {
    let mut cs = vec![0.0f64; jl];

    for i in 0..jl {
        let mut s = 1.0;
        let i_f = (i + 1) as f64;
        for j in 0..jl {
            if i != j {
                let j_f = (j + 1) as f64;
                s = s * (2.0 * j_f - 1.0).powf(2.0)
                    / ((2.0 * i_f - 1.0).powf(2.0) - (2.0 * j_f - 1.0).powf(2.0));
            }
        }
        cs[i] = s.abs() * (-1.0f64).powf(i_f - 1.0) / (2.0 * i_f - 1.0);
    }

    cs
}

// vamx:最大速度
// dt: 采样间隔
// f: 主频
pub fn is_stable(vmax: f64, dt: f64, dx: f64, dz: f64, cs: &[f64]) -> bool {
    let mut d = 0.0;
    for n in 0..cs.len() {
        let n_f = (n + 1) as f64;
        d += cs[n] + (-1.0f64).powf(n_f - 1.0);
    }
    let mut r = vmax.powf(2.0) * dt.powf(2.0) * (1.0 / dx.powf(2.0) + 1.0 / dz.powf(2.0)) * d;
    dbg!(r);
    r <= 1.0 && r >= 0.0
}

fn vp_vs_p(cur_dir_str: &str, para: &io::ModelPara) -> (io::grd::Grd, io::grd::Grd, io::grd::Grd) {
    let vp_grd = format!("{}\\{}", cur_dir_str, para.vp_grd());
    let vs_grd = format!("{}\\{}", cur_dir_str, para.vs_grd());
    let pp_grd = format!("{}\\{}", cur_dir_str, para.pp_grd());

    let vp_grd = io::grd::from_grd_file(&vp_grd);
    let vs_grd = io::grd::from_grd_file(&vs_grd);
    let pp_grd = io::grd::from_grd_file(&pp_grd);

    (vp_grd, vs_grd, pp_grd)
}

// path: PARAMETER.TXT
pub fn run<P: AsRef<Path>>(path: P) {
    let path: &Path = path.as_ref();

    let model_para = match path.is_file() {
        true => io::ModelPara::from_parameter_txt(path),
        false => panic!("wrong PARAMETER.TXT path"),
    };

    let cur_dir = path.parent().unwrap();
    let cur_dir_str = cur_dir.to_str().unwrap();

    //读取参数
    let (vp_grd, vs_grd, pp_grd) = vp_vs_p(cur_dir_str, &model_para);

    let (dx, dz, nx, nz) = (
        *vp_grd.x_size(),
        *vp_grd.y_size(),
        *vp_grd.cols() as usize,
        *vp_grd.rows() as usize,
    );
    let (vpmax, vpmin, vsmax, vsmin) = (
        *vp_grd.z_max(),
        vp_grd.z_min(),
        *vs_grd.z_max(),
        vs_grd.z_min(),
    );
    let (vp, vs, pp) = (vp_grd.data(), vs_grd.data(), pp_grd.data());

    let (sx, sz) = (model_para.source_x(), model_para.source_z());
    let (dt, nt, fm, t0) = (
        *model_para.dt(),
        *model_para.points(),
        *model_para.fm(),
        *model_para.delay(),
    );

    let (jl, pml_h) = (*model_para.diff_order() / 2, *model_para.pml_h());

    let cn = cs(jl);

    if !is_stable(vpmax, dt, dx, dz, &cn) {
        panic!("未达到稳定性条件");
    }

    let nx_pml = nx + 2 * pml_h;
    let nz_pml = nz + 2 * pml_h; //后面修改成面波

    let (lambda, miu) = crate::lame::lame(vp, vs, pp);
    let wave = crate::wavelet::ricker(fm, t0, dt, nt);

    // PML吸收系数
    let r = 1e-6f64;

    let d_pml_x0 = (1.0 / r).ln() * 3.0 * vsmax / (2.0 * nx_pml as f64);
    let mut d_pml_x = Array2::<f64>::zeros((nz_pml, nx_pml));
    let mat1 = Array2::<f64>::ones((nz_pml, 1).f());
    let mat2 = Array2::<f64>::from_shape_fn((1, pml_h).f(), |ix| {
        let j = ix.1;
        ((pml_h - j) as f64 / pml_h as f64).powf(2.0)
    });

    d_pml_x
        .slice_mut(s![.., 0..pml_h])
        .assign(&mat1.dot(&(d_pml_x0 * mat2)));
    let d_pml_x_tmp = d_pml_x.slice(s![.., 0..pml_h; -1]).to_owned(); //解决
    d_pml_x
        .slice_mut(s![.., nx + pml_h..nx_pml])
        .assign(&d_pml_x_tmp);

    let d_pml_z0 = d_pml_x0;
    let mut d_pml_z = Array2::<f64>::zeros((nz_pml, nx_pml));
    let mat1 = Array2::<f64>::ones((1, nx_pml).f());
    let mat2 = Array2::<f64>::from_shape_fn((pml_h, 1).f(), |ix| {
        let j = ix.0;
        ((pml_h - j) as f64 / pml_h as f64).powf(2.0)
    });

    d_pml_z
        .slice_mut(s![0..pml_h, ..])
        .assign(&((d_pml_z0 * mat2).dot(&mat1)));
    let d_pml_z_tmp = d_pml_z.slice(s![0..pml_h;-1, ..]).to_owned(); //解决
    d_pml_z
        .slice_mut(s![nz + pml_h..nz_pml, ..])
        .assign(&d_pml_z_tmp);

    //波场计算
    let lame_c = lambda.clone() + 2.0 * miu.clone();

    let c1 = (2.0 - dt * &d_pml_x).div(2.0 + dt * &d_pml_x);
    let c2 = (2.0 - dt * &d_pml_z).div(2.0 + dt * &d_pml_z);
    let c3 = (2.0 * dt).div(2.0 + dt * &d_pml_x).div(pp).div(dx);
    let c4 = (2.0 * dt).div(2.0 + dt * &d_pml_z).div(pp).div(dz);
    let c5 = (2.0 * dt).div(2.0 + dt * &d_pml_x) * (&lame_c).div(dx);
    let c6 = (2.0 * dt).div(2.0 + dt * &d_pml_z) * (&lambda).div(dz);
    let c7 = (2.0 * dt).div(2.0 + dt * &d_pml_x) * (&lambda).div(dx);
    let c8 = (2.0 * dt).div(2.0 + dt * &d_pml_z) * (&lame_c).div(dz);
    let c9 = (2.0 * dt).div(2.0 + dt * &d_pml_x) * (&miu).div(dx);
    let c0 = (2.0 * dt).div(2.0 + dt * d_pml_z) * (&miu).div(dz);

    let nz_pml_jl = nz_pml + 2 * jl;
    let nx_pml_jl = nx_pml + 2 * jl;

    let z_nodes_pml = jl..(nz_pml_jl - jl);
    let x_nodes_pml = jl..(nx_pml_jl - jl);
    let z_nodes = (jl + pml_h)..(jl + pml_h + nz);
    let x_nodes = (jl + pml_h)..(jl + pml_h + nx);

    let nsrcx = jl + pml_h + sx;
    let nsrcz = jl + pml_h + sz;

    let zeros = Array2::<f64>::zeros((nz_pml_jl, nx_pml_jl));
    let vxt = [zeros.clone(), zeros.clone()];
    let vxx = [zeros.clone(), zeros.clone()];
    let vxz = [zeros.clone(), zeros.clone()];
    let vzt = [zeros.clone(), zeros.clone()];
    let vzx = [zeros.clone(), zeros.clone()];
    let vzz = [zeros.clone(), zeros.clone()];
    let txxt = [zeros.clone(), zeros.clone()];
    let txxx = [zeros.clone(), zeros.clone()];
    let txxz = [zeros.clone(), zeros.clone()];
    let tzzz = [zeros.clone(), zeros.clone()];
    let txzt = [zeros.clone(), zeros.clone()];
    let txzz = [zeros.clone(), zeros];
    let p_sum = std::f64::NAN * Array2::<f64>::ones((nz_pml, nx_pml));
    let p = vec![std::f64::NAN * Array2::<f64>::ones((nz, nx)); nt];
}
