use crate::numerics::numerov::numerov_uniform_grid_renorm;

/// Radial Schrödinger equation for u(r) = r R(r) in units ħ=1.
/// With reduced mass mu:
///   - (1/2mu) u''(r) + [V(r) + l(l+1)/(2mu r^2)] u(r) = E u(r)
///
/// Rearranged to Numerov form:
///   u''(r) = 2mu * (V_eff(r) - E) * u(r)
///
/// Boundary near r=0: u(r) ~ r^{l+1}.
pub fn solve_radial_shooting<V>(
    r_min: f64,
    r_max: f64,
    n: usize,
    mu: f64,
    l: usize,
    e: f64,
    v: &V,
) -> (Vec<f64>, Vec<f64>)
where
    V: Fn(f64) -> f64,
{
    let dr = (r_max - r_min) / (n as f64 - 1.0);

    let lfac = (l as f64) * (l as f64 + 1.0);
    let k = |r: f64| {
        let eps = 1e-12;
        let rr = r.max(eps);
        let v_eff = v(rr) + lfac / (2.0 * mu * rr * rr);
        2.0 * mu * (v_eff - e)
    };

    // initial u(r) ~ r^{l+1}
    let u0 = r_min.powi((l as i32) + 1);
    // second point using same scaling
    let r1 = r_min + dr;
    let u1 = r1.powi((l as i32) + 1);

    let u = numerov_uniform_grid_renorm(r_min, dr, n, u0, u1, k, 200, 1e6);

    let r: Vec<f64> = (0..n).map(|i| r_min + i as f64 * dr).collect();
    (r, u)
}

/// Count nodes (sign changes) in u(r) excluding tiny near-zero noise.
pub fn count_nodes(u: &[f64]) -> usize {
    let mut nodes = 0usize;
    let mut prev = u[0];
    for &val in u.iter().skip(1) {
        if prev == 0.0 {
            prev = val;
            continue;
        }
        if val == 0.0 {
            continue;
        }
        if prev.signum() != val.signum() {
            nodes += 1;
        }
        prev = val;
    }
    nodes
}

/// Residual for shooting: u(r_max).
pub fn residual_at_rmax<V>(
    r_min: f64,
    r_max: f64,
    n: usize,
    mu: f64,
    l: usize,
    e: f64,
    v: &V,
) -> f64
where
    V: Fn(f64) -> f64,
{
    let (_r, u) = solve_radial_shooting(r_min, r_max, n, mu, l, e, v);
    *u.last().unwrap()
}
