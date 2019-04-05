use std::{
    fs,
    fs::File,
    io::{prelude::*, BufReader, Error, SeekFrom},
    mem,
    path::Path,
};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use ndarray::{Array2, ShapeBuilder};

use crate::io::util::*;

use crate::io::grd;

#[derive(Debug, Clone, Getters)]
#[get = "pub"]
pub struct Grd {
    mark: String, //file type
    rows: i32,    //number of rows
    cols: i32,
    xll: f64, //X coordinate of the lower left corner of the grid
    yll: f64,
    x_size: f64, //spacing between adjacent nodes in the X direction (between columns)
    y_size: f64,
    z_min: f64, //minimum Z value within the grid
    z_max: f64,
    pub(crate) data: Array2<f64>, //C order
}

impl Grd {

}

pub fn from_grd_file(filename: &str) -> Grd {
    let mut grd_file = BufReader::new(File::open(filename).expect("error in opening grd file"));

    let mark = grd_file.read_str(4);

    match mark.as_str() {
        "DSAA" => grd::from_ascii_file(grd_file, mark),
        _ => grd::from_binary_file(grd_file, mark),
    }
}

fn from_binary_file(mut grd_file: BufReader<File>, mark: String) -> Grd {
    grd_file
        .seek(SeekFrom::Start(20))
        .expect("Error in seek to 20");

    let rows = grd_file.read_i32::<LittleEndian>().unwrap();

    let cols = grd_file.read_i32::<LittleEndian>().unwrap();

    let xll = grd_file.read_f64::<LittleEndian>().unwrap();

    let yll = grd_file.read_f64::<LittleEndian>().unwrap();

    let x_size = grd_file.read_f64::<LittleEndian>().unwrap();

    let y_size = grd_file.read_f64::<LittleEndian>().unwrap();

    let z_min = grd_file.read_f64::<LittleEndian>().unwrap();

    let z_max = grd_file.read_f64::<LittleEndian>().unwrap();

    grd_file
        .seek(SeekFrom::Start(100))
        .expect("Error in seek file to 100");

    let nrows = rows as usize;
    let ncols = cols as usize;

    let mut data_vec = vec![0.0; nrows * ncols];

    for i in 0..nrows {
        for j in 0..ncols {
            data_vec[i * ncols + j] = grd_file.read_f64::<LittleEndian>().unwrap();
        }
    }

    let data_r = Array2::from_shape_vec((nrows, ncols).f(), data_vec).unwrap();
    let data = Array2::from_shape_fn((nrows, ncols).f(), |ix| {
        let i = ix.0;
        let j = ix.1;
        data_r[[nrows - i - 1, j]]
    });

    Grd {
        mark,
        rows,
        cols,
        xll,
        yll,
        x_size,
        y_size,
        z_min,
        z_max,
        data,
    }
}

fn from_ascii_file(grd_file: BufReader<File>, _mark: String) -> Grd {
    let buf_reader = BufReader::new(grd_file);

    let mut lines = buf_reader.lines();

    let mark = grd::process_error(lines.next());
    let (cols, rows) = grd::process_split::<i32>(&grd::process_error(lines.next()));
    let (xll, xend) = grd::process_split::<f64>(&grd::process_error(lines.next()));
    let (yll, yend) = grd::process_split::<f64>(&grd::process_error(lines.next()));
    let (z_min, z_max) = grd::process_split::<f64>(&grd::process_error(lines.next()));
    let x_size = (xend - xll) / (cols - 1) as f64;
    let y_size = (yend - yll) / (rows - 1) as f64;

    let mut r: Vec<_> = lines.map(|f| f.unwrap()).collect();
    r.reverse();
    let data_str = r.join(" ");
    let data_intern: Vec<f64> = data_str
        .split_whitespace()
        .map(|x| x.parse::<f64>().unwrap())
        .collect();

    let data = Array2::from_shape_vec((rows as usize, cols as usize).f(), data_intern).unwrap(); // row-major matrix.

    Grd {
        mark,
        rows,
        cols,
        xll,
        yll,
        x_size,
        y_size,
        z_min,
        z_max,
        data,
    }
}

fn process_error(line: Option<Result<String, Error>>) -> String {
    line.unwrap().unwrap()
}

fn process_split<U>(s: &str) -> (U, U)
where
    U: std::str::FromStr,
    <U as std::str::FromStr>::Err: std::fmt::Debug,
{
    let ss: Vec<&str> = s
        .trim()
        .split(char::is_whitespace)
        .filter(|item| *item != "")
        .collect();
    let v1 = ss[0].parse::<U>().unwrap();
    let v2 = ss[1].parse::<U>().unwrap();
    (v1, v2)
}
