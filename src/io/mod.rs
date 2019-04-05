pub mod grd;
pub mod util;

use std::{
    path::Path,
    fs,
    io::{prelude::*, BufReader},
};

// PARAMETER.TXT
#[derive(Debug, Clone, Getters, MutGetters)]
#[get = "pub"]
pub struct ModelPara {
    diff_order: usize,  // 差分阶数
    pml_h: usize, // PML边界网格点数
    dt: f64, // 采样间隔
    points: usize, // 采样点数
    fm: f64, // 震源主频
    delay: f64, // 延迟时(ms)
    source_x: usize, // 震源X方向节点序号
    source_z: usize, // Z方向
    receiver_z: usize, // 检波器Z方向
    source_type: usize, //1胀缩震源，2垂直震源，3水平震源
    vp_grd: String,
    vs_grd: String,
    pp_grd: String,
    wavelet_bln: String,
    cdp_x2: String,
    cdp_z2: String,
    wave_field_x: String,
    wave_field_z: String,
    self_start: usize,
    self_end: usize,
}

impl ModelPara {
    fn to_str_vec(&self) -> Vec<String> {
        let mut str_vec = Vec::new();
        str_vec.push(format!(
            "{},{}",
            self.diff_order, self.pml_h
        ));
        str_vec.push(format!(
            "{},{},{},{}",
            self.dt, self.points, self.fm, self.delay
        ));
        str_vec.push(format!(
            "{},{},{}",
            self.source_x, self.source_z, self.receiver_z
        ));
        str_vec.push(format!("{}", self.source_type));
        str_vec.push(format!("{}", self.vp_grd));
        str_vec.push(format!("{}", self.vs_grd));
        str_vec.push(format!("{}", self.pp_grd));
        str_vec.push(format!("{}", self.wavelet_bln));
        str_vec.push(format!("{}", self.cdp_x2));
        str_vec.push(format!("{}", self.cdp_z2));
        str_vec.push(format!("{}", self.wave_field_x));
        str_vec.push(format!("{}", self.wave_field_z));
        str_vec.push(format!("{},{}", self.self_start, self.self_end));
        str_vec
    }

    pub fn from_parameter_txt<P: AsRef<Path>>(path: P) -> ModelPara {
        let file = BufReader::new(fs::File::open(path).expect("error in open parameter.txt"));

        let mut lines: Vec<_> = file.lines().map(|l| l.unwrap()).collect();
        let v1: Vec<_> = lines[1]
            .split(',')
            .map(|s| s.parse::<usize>().unwrap())
            .collect();
        let (diff_order, pml_h) = (v1[0], v1[1]);

        let v3: Vec<_> = lines[3]
            .split(',')
            .map(|s| s.parse::<f64>().unwrap())
            .collect();
        let (dt, points, fm, delay) = (v3[0], v3[1] as usize, v3[2], v3[3]);

        let v5: Vec<_> = lines[5]
            .split(',')
            .map(|s| s.parse::<usize>().unwrap())
            .collect();
        let (source_x, source_z, receiver_z) = (v5[0], v5[1], v5[2]);

        let v7: Vec<_> = lines[7]
            .split(',')
            .map(|s| s.parse::<usize>().unwrap())
            .collect();
        let source_type = v7[0];
        // let save_or_not = 1;
        let vp_grd = lines[9].clone();
        let vs_grd = lines[11].clone();
        let pp_grd = lines[13].clone();
        let wavelet_bln = lines[15].clone();
        let cdp_x2 = lines[17].clone();
        let cdp_z2 = lines[19].clone();
        let wave_field_x = lines[21].clone();
        let wave_field_z = lines[23].clone();
        let v25: Vec<_> = lines[25]
            .split(',')
            .map(|s| s.parse::<usize>().unwrap())
            .collect();
        let (self_start, self_end) = (v25[0], v25[1]);

        ModelPara {
            diff_order,
            pml_h,
            dt,
            points,
            fm,
            delay,
            source_x,
            source_z,
            receiver_z,
            source_type,
            vp_grd,
            vs_grd,
            pp_grd,
            wavelet_bln,
            cdp_x2,
            cdp_z2,
            wave_field_x,
            wave_field_z,
            self_start,
            self_end,
        }
    }

    pub fn write(&self, dir_name: &str) {
        let path_str = format!("{}\\PARAMETER.txt", dir_name);
        let path = Path::new(&path_str);
        let mut file = fs::File::create(path).expect("error in create parameter file");
        let model_template = vec![
            "!差分阶数，吸收边界厚度",
            "!采样间隔，采样点数，震源子波主频，震源子波延迟时",
            "!震源所在X、Z方向节点序号，检波点所在Z方向节点序号",
            "!点震源的类型：1胀缩震源，2垂直震源，3水平震源",
            "!模型的纵波速度",
            "!模型的横波速度",
            "!模型的密度",
            "!保存震源子波的文件名",
            "!保存共炮点记录水平分量的文件名",
            "!保存共炮点记录垂直分量的文件名",
            "!保存波场水平分量的文件名",
            "!保存波场垂直分量的文件名",
            "!自激自收区域的起点、终点位置",
        ];
        // let model_template: Vec<_> = model_template.iter().map(|x| x.to_string()).collect();
        // dbg!(&model_template);
        let model_values = self.to_str_vec();
        assert!(model_template.len() == model_values.len());

        for i in 0..model_template.len() {
            write!(file, "{}\r\n", model_template[i]);
            write!(file, "{}\r\n", model_values[i]);
        }
    }

    pub fn new_with_default(model_prefix: &str) -> Self {
        ModelPara {
            diff_order: 10,
            pml_h: 100,
            dt: 0.0002,
            points: 5001,
            fm: 20.0,
            delay: 0.05,
            source_x: 1,
            source_z: 1,
            receiver_z: 1,
            source_type: 2,
            vp_grd: format!("{}vp.grd", model_prefix),
            vs_grd: format!("{}vs.grd", model_prefix),
            pp_grd: format!("{}vp.grd", model_prefix),
            wavelet_bln: "wavelet.bln".to_string(),
            cdp_x2: format!("{}-point-x2.cdp", model_prefix),
            cdp_z2: format!("{}-point-z2.cdp", model_prefix),
            wave_field_x: format!("{}-point-X.dat", model_prefix), // 波场快照
            wave_field_z: format!("{}-point-Z.dat", model_prefix),
            self_start: 1,
            self_end: 21,
        }
    }
}