use anyhow::{Context, Result};
use clap::Parser;
use plotters::prelude::*;
use std::fs::create_dir_all;

use rust4physics::quantum::wavepacket::free_gaussian_density;

#[derive(Parser, Debug)]
struct Args {
    /// Output directory for frames
    #[arg(long, default_value = "out/wp_frames")]
    outdir: String,

    /// Number of frames
    #[arg(long, default_value_t = 80)]
    frames: usize,

    /// Grid min x
    #[arg(long, default_value_t = -40.0)]
    xmin: f64,

    /// Grid max x
    #[arg(long, default_value_t = 40.0)]
    xmax: f64,

    /// Number of grid points
    #[arg(long, default_value_t = 2000)]
    n: usize,

    /// Total simulated time
    #[arg(long, default_value_t = 12.0)]
    tmax: f64,

    /// Initial packet center x0
    #[arg(long, default_value_t = -12.0)]
    x0: f64,

    /// Initial width sigma0
    #[arg(long, default_value_t = 1.2)]
    sigma0: f64,

    /// Central wavenumber k0 (group velocity v = k0 for ħ=m=1)
    #[arg(long, default_value_t = 2.0)]
    k0: f64,
}

fn main() -> Result<()> {
    let args = Args::parse();
    create_dir_all(&args.outdir).context("create outdir")?;

    let dx = (args.xmax - args.xmin) / (args.n as f64 - 1.0);

    for frame in 0..args.frames {
        let t = args.tmax * (frame as f64) / ((args.frames - 1).max(1) as f64);

        let mut y = Vec::with_capacity(args.n);
        let mut ymax: f64 = 0.0;

        for i in 0..args.n {
            let x = args.xmin + i as f64 * dx;
            let rho = free_gaussian_density(x, t, args.x0, args.sigma0, args.k0);
            ymax = ymax.max(rho);
            y.push((x, rho));
        }

        let path = format!("{}/frame_{:04}.png", args.outdir, frame);
        plot_frame(&path, &y, args.xmin, args.xmax, ymax, t, args.k0)?;
        println!("Wrote {}", path);
    }

    Ok(())
}

fn plot_frame(path: &str, data: &[(f64, f64)], xmin: f64, xmax: f64, ymax: f64, t: f64, k0: f64) -> Result<()> {
    let root = BitMapBackend::new(path, (1100, 520)).into_drawing_area();
    root.fill(&WHITE)?;

    let title = format!("Free Gaussian wavepacket: |ψ(x,t)|²   t={:.3}   (ħ=1,m=1, k0={:.2})", t, k0);

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(title, ("sans-serif", 22))
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(xmin..xmax, 0.0..(ymax * 1.15))?;

    chart.configure_mesh().x_desc("x").y_desc("|ψ|²").draw()?;
    chart.draw_series(LineSeries::new(data.iter().copied(), &BLUE))?;

    root.present()?;
    Ok(())
}
