use anyhow::{Context, Result};
use clap::Parser;
use plotters::prelude::*;
use std::fs::{create_dir_all, File};
use std::io::Write;

use rust4physics::numerics::root::bisect;
use rust4physics::quantum::medium::{bjorken_temperature, mu_debye};
use rust4physics::quantum::potentials::screened_cornell;
use rust4physics::quantum::radial::residual_at_rmax;

#[derive(Parser, Debug)]
struct Args {
    #[arg(long, default_value = "out/bjorken")]
    outdir: String,

    /// Reduced mass μ (GeV)
    #[arg(long, default_value_t = 2.4)]
    mu: f64,

    /// Cornell alpha, sigma, C
    #[arg(long, default_value_t = 0.52)]
    alpha: f64,
    #[arg(long, default_value_t = 0.18)]
    sigma: f64,
    #[arg(long, default_value_t = 0.0)]
    c: f64,

    /// Debye coefficient mu(T)=a*T
    #[arg(long, default_value_t = 2.5)]
    a_debye: f64,

    /// Bjorken params (GeV units for T)
    #[arg(long, default_value_t = 0.6)]
    tau0: f64,
    #[arg(long, default_value_t = 0.5)]
    t0: f64,
    #[arg(long, default_value_t = 6.0)]
    taumax: f64,

    /// stop at Tc
    #[arg(long, default_value_t = 0.170)]
    tc: f64,

    /// radial grid
    #[arg(long, default_value_t = 1e-4)]
    rmin: f64,
    #[arg(long, default_value_t = 10.0)]
    rmax: f64,
    #[arg(long, default_value_t = 6000)]
    n: usize,

    /// angular momentum l
    #[arg(long, default_value_t = 0)]
    l: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    create_dir_all(&args.outdir).context("create outdir")?;

    // sample tau values
    let n_tau = 80usize;
    let mut rows: Vec<(f64, f64, f64)> = Vec::new(); // (tau, T, E0)

    for i in 0..n_tau {
        let tau = args.tau0 + (args.taumax - args.tau0) * (i as f64) / ((n_tau - 1) as f64);
        let t = bjorken_temperature(tau, args.tau0, args.t0);
        if t < args.tc {
            break;
        }
        let mu_screen = mu_debye(t, args.a_debye);

        let v = |r: f64| screened_cornell(r, args.alpha, args.sigma, mu_screen, args.c);

        // bracket ground-ish state energy in a toy interval
        let mut lo = -2.0;
        let hi =  0.0;

        // ensure bracket sign change
        let mut flo = residual_at_rmax(args.rmin, args.rmax, args.n, args.mu, args.l, lo, &v);
        let fhi = residual_at_rmax(args.rmin, args.rmax, args.n, args.mu, args.l, hi, &v);

        // expand if needed (limited)
        for _ in 0..20 {
            if flo.signum() != fhi.signum() {
                break;
            }
            lo -= 0.5;
            flo = residual_at_rmax(args.rmin, args.rmax, args.n, args.mu, args.l, lo, &v);
        }
        if flo.signum() == fhi.signum() {
            // couldn't find bound root in this toy scan
            rows.push((tau, t, f64::NAN));
            continue;
        }

        let e0 = bisect(lo, hi, |en| residual_at_rmax(args.rmin, args.rmax, args.n, args.mu, args.l, en, &v), 1e-8, 250)
            .map_err(|s| anyhow::anyhow!(s))?;

        rows.push((tau, t, e0));
    }

    // write csv
    let csv = format!("{}/scan.csv", args.outdir);
    let mut f = File::create(&csv).context("write csv")?;
    writeln!(f, "tau,T,E0")?;
    for (tau, t, e) in &rows {
        writeln!(f, "{},{},{}", tau, t, e)?;
    }
    println!("Wrote {}", csv);

    // plot E0 vs tau
    let png = format!("{}/E0_vs_tau.png", args.outdir);
    plot_scan(&png, &rows)?;
    println!("Wrote {}", png);

    Ok(())
}

fn plot_scan(path: &str, rows: &[(f64, f64, f64)]) -> Result<()> {
    let root = BitMapBackend::new(path, (1050, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let xmin = rows.first().map(|r| r.0).unwrap_or(0.0);
    let xmax = rows.last().map(|r| r.0).unwrap_or(1.0);

    // y-range from finite values
    let mut ymin: f64 =  1e9;
    let mut ymax: f64 = -1e9;
    for &(_, _, e) in rows {
        if e.is_finite() {
            ymin = ymin.min(e);
            ymax = ymax.max(e);
        }
    }
    if ymin > ymax {
        ymin = -2.0; ymax = 0.0;
    }
    let pad: f64 = 0.1 * (ymax - ymin).abs().max(1e-3);
    ymin -= pad; ymax += pad;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption("Toy screened Cornell: ground energy vs Bjorken proper time", ("sans-serif", 22))
        .x_label_area_size(45)
        .y_label_area_size(70)
        .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

    chart.configure_mesh().x_desc("tau").y_desc("E0 [GeV]").draw()?;

    chart.draw_series(LineSeries::new(
        rows.iter().filter(|(_,_,e)| e.is_finite()).map(|(tau,_,e)| (*tau, *e)),
        &BLUE,
    ))?;

    root.present()?;
    Ok(())
}
