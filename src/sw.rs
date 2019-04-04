use ndarray::{Array2, Array, ShapeBuilder};

pub fn cs(jl: usize) -> Vec<f64>{
    let mut cs = vec![0.0f64; jl];

    for i in 0..jl {
        let mut s = 1.0;
        let i_f = (i + 1) as f64;
        for j in 0..jl {
            if i != j {
                let j_f = (j + 1) as f64;
                s = s * (2.0 * j_f - 1.0).powf(2.0) / ((2.0*i_f-1.0).powf(2.0) - (2.0*j_f-1.0).powf(2.0));
            }
        }
        cs[i] = s.abs() * (-1.0f64).powf(i_f-1.0)/(2.0*i_f-1.0);
    }

    cs
}

// vamx:最大速度
// dt: 采样间隔
// f: 主频
pub fn is_stable(vmax: f64, dt: f64, dx: f64, dz: f64, cs: &Vec<f64>) -> bool {
    let mut d = 0.0;
    for n in 0..cs.len() {
        let n_f = (n + 1) as f64;
        d += cs[n] + (-1.0f64).powf(n_f-1.0);
    }
    let mut r = vmax.powf(2.0) * dt.powf(2.0) * (1.0 / dx.powf(2.0) + 1.0 / dz.powf(2.0)) * d;
    dbg!(r);
    r <= 1.0 && r >= 0.0 
}


pub fn run(dir: &str) {
    
}