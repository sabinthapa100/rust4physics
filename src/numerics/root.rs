/// Bisection root finder on [lo, hi], requiring f(lo) and f(hi) to have opposite signs.
pub fn bisect<F>(
    mut lo: f64,
    mut hi: f64,
    f: F,
    tol: f64,
    max_iter: usize,
) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let mut flo = f(lo);
    let fhi = f(hi);

    if flo == 0.0 {
        return Ok(lo);
    }
    if fhi == 0.0 {
        return Ok(hi);
    }
    if flo.signum() == fhi.signum() {
        return Err("bisect: need bracket with opposite signs".into());
    }

    for _ in 0..max_iter {
        let mid = 0.5 * (lo + hi);
        let fmid = f(mid);

        if fmid.abs() < tol || (hi - lo).abs() < tol {
            return Ok(mid);
        }

        if fmid.signum() == flo.signum() {
            lo = mid;
            flo = fmid;
        } else {
            hi = mid;
            // fhi not needed; hi update is enough for bisection
        }
    }

    Ok(0.5 * (lo + hi))
}
