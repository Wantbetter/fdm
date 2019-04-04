
pub fn ricker(f0: f64, t0: f64, dt: f64, nt: usize) -> Vec<f64> {
    let mut wavelet = vec![0.0f64; nt];
    let pi = std::f64::consts::PI;
    let e = std::f64::consts::E;
    for i in 0..wavelet.len() {
        let t = dt * i as f64 - t0;
        let combine = (pi * f0 * t).powf(2.0);
        wavelet[i] = (1.0 - 2.0 * combine) * e.powf(-combine);
    } 
    wavelet
}