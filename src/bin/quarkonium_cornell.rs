use anyhow::{Context, Result};
use clap::Parser;
use plotters::prelude::*;
use std::fs::{create_dir_all, File};
use std::io::Write;

use rust4physics::numerics::root::bisect;
use rust4physics::quantum::potentials::cornell;
use rust4physics::quantum::radial::{residual_at_rmax, solve_radial_shooting};

#[derive(Parser, Debug)]
struct Args {
    #[arg(long, default_value = "out/cornell")]
    outdir: String,

    /// Reduced mass μ (GeV) in natural units ħ=c=1 (toy default for bottomonium)
    #[arg(long, default_value_t = 2.4)]
    mu: f64,

    /// Cornell alpha
    #[arg(long, default_value_t = 0.52)]
    alpha: f64,

    /// Cornell sigma (GeV^2) if r in GeV^-1
    #[arg(long, default_value_t = 0.18)]
    sigma: f64,

    /// Constant shift C (GeV)
    #[arg(long, default_value_t = 0.0)]
    c: f64,

    /// r_min
    #[arg(long, default_value_t = 1e-4)]
    rmin: f64,

    /// r_max
    #[arg(long, default_value_t = 10.0)]
    rmax: f64,

    /// grid points
    #[arg(long, default_value_t = 6000)]
    n: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    create_dir_all(&args.outdir).context("create outdir")?;

    let v = |r: f64| cornell(r, args.alpha, args.sigma, args.c);

    // We locate a few low-lying states by scanning for sign changes in residual at rmax.
    // This is a practical “exploration tool”. For production, you’d add node-count control and better bracketing.
    let mut levels: Vec<(usize, usize, f64)> = Vec::new(); // (l, idx, E)

    for &l in &[0usize, 1usize] { // S and P
        let mut brackets = Vec::new();
        let e_min = -2.0;  // GeV (toy)
        let e_max =  0.0;
        let n_scan = 600;

        let mut prev_e = e_min;
        let mut prev_r = residual_at_rmax(args.rmin, args.rmax, args.n, args.mu, l, prev_e, &v);

        for i in 1..=n_scan {
            let e = e_min + (e_max - e_min) * (i as f64) / (n_scan as f64);
            let r = residual_at_rmax(args.rmin, args.rmax, args.n, args.mu, l, e, &v);
            if prev_r.signum() != r.signum() {
                brackets.push((prev_e, e));
            }
            prev_e = e;
            prev_r = r;
        }

        // refine first few
        for (idx, (lo, hi)) in brackets.into_iter().take(4).enumerate() {
            let e = bisect(lo, hi, |en| residual_at_rmax(args.rmin, args.rmax, args.n, args.mu, l, en, &v), 1e-8, 250)
                .map_err(|s| anyhow::anyhow!(s))?;
            levels.push((l, idx + 1, e));

            // plot wavefunction u(r) for this level
            let (rgrid, u) = solve_radial_shooting(args.rmin, args.rmax, args.n, args.mu, l, e, &v);
            let png = format!("{}/state_l{}_k{}.png", args.outdir, l, idx + 1);
            plot_u(&png, &rgrid, &u, l, e)?;
            println!("Wrote {}", png);
        }
    }

    // write levels
    let path = format!("{}/levels.txt", args.outdir);
    let mut f = File::create(&path).context("write levels")?;
    writeln!(f, "Cornell bound-state scan (toy).")?;
    writeln!(f, "mu={} GeV, alpha={}, sigma={} GeV^2, C={} GeV", args.mu, args.alpha, args.sigma, args.c)?;
    writeln!(f, "Columns: l, index(k), E [GeV]")?;
    for (l, k, e) in levels {
        writeln!(f, "{} {} {:.8}", l, k, e)?;
    }
    println!("Wrote {}", path);

    Ok(())
}

fn plot_u(path: &str, r: &[f64], u: &[f64], l: usize, e: f64) -> Result<()> {
    let root = BitMapBackend::new(path, (1050, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let rmin = *r.first().unwrap();
    let rmax = *r.last().unwrap();

    // crude y-range
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for &v in u {
        ymin = ymin.min(v);
        ymax = ymax.max(v);
    }
    let pad = 0.1 * (ymax - ymin).max(1e-6);
    ymin -= pad; ymax += pad;

    let title = format!("Radial u(r) for Cornell (l={}, E={:.4} GeV)", l, e);

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(title, ("sans-serif", 22))
        .x_label_area_size(45)
        .y_label_area_size(70)
        .build_cartesian_2d(rmin..rmax, ymin..ymax)?;

    chart.configure_mesh().x_desc("r").y_desc("u(r)").draw()?;
    chart.draw_series(LineSeries::new(r.iter().copied().zip(u.iter().copied()), &BLUE))?;

    root.present()?;
    Ok(())
}
