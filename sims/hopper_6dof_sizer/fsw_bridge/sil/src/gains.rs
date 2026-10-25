//! Loads the LQR gain schedule written by `export_gains.m`.

use std::{fs, io, path::Path};

use nalgebra::{Const, Dyn, MatrixXx1, MatrixXx4, OMatrix};

use crate::control::lqr::LqrController;

/// Reads a headerless CSV of numbers into (rows, row-major data), checking
/// every row has `cols` values.
fn read_csv(path: &Path, cols: usize) -> io::Result<(usize, Vec<f64>)> {
    let text = fs::read_to_string(path)
        .map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    let mut data = Vec::new();
    let mut rows = 0;
    for (i, line) in text.lines().filter(|l| !l.trim().is_empty()).enumerate() {
        let row: Vec<f64> = line
            .split(',')
            .map(|v| v.trim().parse::<f64>())
            .collect::<Result<_, _>>()
            .map_err(|e| invalid(format!("{} line {}: {e}", path.display(), i + 1)))?;
        if row.len() != cols {
            return Err(invalid(format!(
                "{} line {}: expected {cols} values, got {}",
                path.display(),
                i + 1,
                row.len()
            )));
        }
        data.extend(row);
        rows += 1;
    }
    Ok((rows, data))
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

fn load<const C: usize>(dir: &Path, name: &str) -> io::Result<OMatrix<f64, Dyn, Const<C>>> {
    let (rows, data) = read_csv(&dir.join(format!("{name}.csv")), C)?;
    Ok(OMatrix::<f64, Dyn, Const<C>>::from_row_slice_generic(Dyn(rows), Const::<C>, &data))
}

/// Builds luna's `LqrController` from the exported tables.
///
/// `negate_k2` passes `-K2grid`. luna's `lqr.rs` computes
/// `u = unom - K1*dx + K2` while the Simulink controller computes
/// `u = unom - K1*dx - K2`, so negating the table reproduces Simulink.
pub fn load_lqr(dir: &Path, negate_k2: bool) -> io::Result<LqrController> {
    let tgrid: MatrixXx1<f64> = load::<1>(dir, "tgrid")?;
    let xref = load::<13>(dir, "xref")?;
    let k1flat = load::<52>(dir, "k1flat")?;
    let mut k2grid: MatrixXx4<f64> = load::<4>(dir, "k2grid")?;
    let unom: MatrixXx4<f64> = load::<4>(dir, "unom")?;
    if negate_k2 {
        k2grid = -k2grid;
    }
    Ok(LqrController::new(tgrid, xref, k1flat, k2grid, unom))
}
