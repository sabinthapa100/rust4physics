use anyhow::{Context, Result};
use clap::Parser;
use plotters::prelude::*;
use std::fs::{create_dir_all, File};
use std::io::Write;

use rust4physics::numerics::{numerov::numerov_uniform_grid_renorm, root::bisect};
use rust4physics::quantum::potentials::harmonic_oscillator;

#[derive(Parser, Debug)]
struct Args {
    /// Output directory for data + plots
    #[arg(long, default_value = "out/ho")]
    outdir: String,

    /// Domain max (we integrate from 0..xmax, using parity)
    #[arg(long, default_value_t = 8.0)]
    xmax: f64,

    /// Number of grid points on [0, xmax]
    #[arg(long, default_value_t = 4000)]
    n: usize,

    /// HO omega
    #[arg(long, default_value_t = 1.0)]
    omega: f64,
}

/// In ħ=1, m=1 units: ψ''(x) = 2 (V(x) - E) ψ(x)
fn k_tise(x: f64, e: f64, omega: f64) -> f64 {
    2.0 * (e - harmonic_oscillator(x, omega))
}

/// Even-parity shooting from x=0 to x=xmax.
/// Uses Taylor expansion to seed ψ(dx) given ψ(0)=1, ψ'(0)=0.
/// Residual is ψ(xmax) (sign changes bracket eigenvalues).
fn residual_even(e: f64, dx: f64, n: usize, omega: f64) -> f64 {
    let k = |x: f64| k_tise(x, e, omega);

    let psi0 = 1.0;
    // ψ''(0)=2(V(0)-E)ψ(0)=-2E
    let psi1 = psi0 + 0.5 * (-2.0 * e) * dx * dx;

    let psi = numerov_uniform_grid_renorm(0.0, dx, n, psi0, psi1, k, 200, 1e6);
    *psi.last().unwrap()
}

fn main() -> Result<()> {
    let args = Args::parse();
    create_dir_all(&args.outdir).context("create outdir")?;

    let dx = args.xmax / (args.n as f64 - 1.0);

    // Find first few even-parity eigenvalues: n=0,2,4,... in full HO labeling.
    // In HO: E_n = ω (n + 1/2). Even parity corresponds to n even.
    let targets = [0.5, 2.5, 4.5, 6.5].map(|x| x * args.omega);
    let mut energies = Vec::new();

    for &guess in &targets {
        let mut lo = guess - 0.4 * args.omega;
        let mut hi = guess + 0.4 * args.omega;

        let mut flo = residual_even(lo, dx, args.n, args.omega);
        let mut fhi = residual_even(hi, dx, args.n, args.omega);

        for _ in 0..50 {
            if flo.signum() != fhi.signum() {
                break;
            }
            lo -= 0.2 * args.omega;
            hi += 0.2 * args.omega;
            flo = residual_even(lo, dx, args.n, args.omega);
            fhi = residual_even(hi, dx, args.n, args.omega);
        }
        if flo.signum() == fhi.signum() {
            anyhow::bail!("Failed to bracket eigenvalue near guess={}", guess);
        }

        let e = bisect(lo, hi, |en| residual_even(en, dx, args.n, args.omega), 1e-11, 300)
            .map_err(|s| anyhow::anyhow!(s))?;
        energies.push(e);
    }

    // Save energies
    let energies_path = format!("{}/energies.txt", args.outdir);
    let mut f = File::create(&energies_path).context("write energies")?;
    writeln!(f, "Harmonic oscillator (ħ=1,m=1) omega={}", args.omega)?;
    writeln!(f, "Even parity states only (n=0,2,4,6 in full HO labeling)")?;
    for (i, e) in energies.iter().enumerate() {
        let n_full = 2 * i;
        let expected = args.omega * (n_full as f64 + 0.5);
        writeln!(f, "n_full={}  E={:.12}  expected={:.12}", n_full, e, expected)?;
    }
    println!("Wrote {}", energies_path);

    // Build ground state ψ(x) on [-xmax, xmax] by even symmetry
    let e0 = energies[0];
    let k = |x: f64| k_tise(x, e0, args.omega);

    let psi0 = 1.0;
    let psi1 = psi0 + 0.5 * (-2.0 * e0) * dx * dx;
    let psi_pos = numerov_uniform_grid_renorm(0.0, dx, args.n, psi0, psi1, k, 200, 1e6);

    // Mirror to negative x
    let mut x_full = Vec::with_capacity(2 * args.n - 1);
    let mut psi_full = Vec::with_capacity(2 * args.n - 1);

    for i in (1..args.n).rev() {
        x_full.push(-(i as f64) * dx);
        psi_full.push(psi_pos[i]);
    }
    for i in 0..args.n {
        x_full.push(i as f64 * dx);
        psi_full.push(psi_pos[i]);
    }

    // Normalize ψ on full domain with trapezoid
    let mut norm = 0.0;
    for i in 0..(psi_full.len() - 1) {
        let a = psi_full[i] * psi_full[i];
        let b = psi_full[i + 1] * psi_full[i + 1];
        norm += 0.5 * (a + b) * dx;
    }
    let norm_sqrt = norm.sqrt();
    for v in psi_full.iter_mut() {
        *v /= norm_sqrt;
    }

    // Save CSV
    let csv_path = format!("{}/ground.csv", args.outdir);
    let mut g = File::create(&csv_path).context("write csv")?;
    writeln!(g, "x,psi,V")?;
    for (&x, &psi) in x_full.iter().zip(psi_full.iter()) {
        let v = harmonic_oscillator(x, args.omega);
        writeln!(g, "{},{},{}", x, psi, v)?;
    }
    println!("Wrote {}", csv_path);

    // Plot
    let png_path = format!("{}/ground.png", args.outdir);
    plot_ground(&png_path, &x_full, &psi_full, args.omega)?;
    println!("Wrote {}", png_path);

    Ok(())
}

fn plot_ground(path: &str, x: &[f64], psi: &[f64], omega: f64) -> Result<()> {
    let xmin = *x.first().unwrap();
    let xmax = *x.last().unwrap();

    let root = BitMapBackend::new(path, (1100, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    let (top, bottom) = root.split_vertically(410);

    // ψ
    {
        let mut chart = ChartBuilder::on(&top)
            .margin(20)
            .caption("Harmonic Oscillator ground state ψ(x) (Numerov, ħ=1,m=1)", ("sans-serif", 22))
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(xmin..xmax, -0.8..0.8)?;

        chart.configure_mesh().x_desc("x").y_desc("ψ(x)").draw()?;
        chart.draw_series(LineSeries::new(
            x.iter().copied().zip(psi.iter().copied()),
            &BLUE,
        ))?;
    }

    // V
    {
        let mut chart = ChartBuilder::on(&bottom)
            .margin(20)
            .caption("Potential V(x) = 0.5 ω² x²", ("sans-serif", 18))
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(xmin..xmax, 0.0..(0.5 * omega * omega * xmax * xmax).min(30.0))?;

        chart.configure_mesh().x_desc("x").y_desc("V(x)").draw()?;
        chart.draw_series(LineSeries::new(
            x.iter().map(|&xx| (xx, harmonic_oscillator(xx, omega))),
            &RED,
        ))?;
    }

    root.present()?;
    Ok(())
}
