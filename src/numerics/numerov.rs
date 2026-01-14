/// Numerov integrates y''(x) = k(x) y(x) on a uniform grid.
/// Returns y values on the grid, given initial y[0] and y[1].
pub fn numerov_uniform_grid<K>(
    x0: f64,
    dx: f64,
    n: usize,
    y0: f64,
    y1: f64,
    k: K,
) -> Vec<f64>
where
    K: Fn(f64) -> f64,
{
    assert!(n >= 2, "need at least 2 grid points");
    let mut y = vec![0.0_f64; n];
    y[0] = y0;
    y[1] = y1;

    let h2 = dx * dx;

    for i in 1..(n - 1) {
        let xm1 = x0 + (i as f64 - 1.0) * dx;
        let xi = x0 + (i as f64) * dx;
        let xp1 = x0 + (i as f64 + 1.0) * dx;

        let km1 = k(xm1);
        let ki = k(xi);
        let kp1 = k(xp1);

        let a = 1.0 + h2 * kp1 / 12.0;
        let b = 2.0 * (1.0 - 5.0 * h2 * ki / 12.0) * y[i];
        let c = (1.0 + h2 * km1 / 12.0) * y[i - 1];

        y[i + 1] = (b - c) / a;
    }

    y
}

/// Same as Numerov, but with periodic renormalization to prevent overflow.
/// This is helpful for bound-state shooting, where a tiny admixture of the growing solution can blow up.
pub fn numerov_uniform_grid_renorm<K>(
    x0: f64,
    dx: f64,
    n: usize,
    y0: f64,
    y1: f64,
    k: K,
    renorm_every: usize,
    renorm_threshold: f64,
) -> Vec<f64>
where
    K: Fn(f64) -> f64,
{
    assert!(n >= 2);
    let mut y = vec![0.0_f64; n];
    y[0] = y0;
    y[1] = y1;

    let h2 = dx * dx;

    for i in 1..(n - 1) {
        let xm1 = x0 + (i as f64 - 1.0) * dx;
        let xi = x0 + (i as f64) * dx;
        let xp1 = x0 + (i as f64 + 1.0) * dx;

        let km1 = k(xm1);
        let ki = k(xi);
        let kp1 = k(xp1);

        let a = 1.0 + h2 * kp1 / 12.0;
        let b = 2.0 * (1.0 - 5.0 * h2 * ki / 12.0) * y[i];
        let c = (1.0 + h2 * km1 / 12.0) * y[i - 1];

        y[i + 1] = (b - c) / a;

        if renorm_every > 0 && i % renorm_every == 0 {
            let scale = y[i].abs().max(y[i + 1].abs());
            if scale > renorm_threshold {
                y[i] /= scale;
                y[i + 1] /= scale;
            }
        }
    }

    y
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn numerov_matches_sin_for_constant_negative_k() {
        // y'' + w^2 y = 0  => y = sin(w x) satisfies y(0)=0, y'(0)=w
        let w: f64 = 2.0;
        let k = |_x: f64| w * w;
        
        let x0 = 0.0;
        let dx: f64 = 1e-3;
        let n = 10_000;

        // y(0)=0, y(dx)=sin(w dx)
        let y0 = 0.0;
        let y1 = (w * dx).sin();

        let y = numerov_uniform_grid(x0, dx, n, y0, y1, k);

        let idx = 9000usize;
        let x = x0 + idx as f64 * dx;

        let exact = (w * x).sin();
        assert_relative_eq!(y[idx], exact, max_relative = 5e-4);
    }
}
