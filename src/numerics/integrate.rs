/// Composite Simpson integration on a uniform grid.
/// `y` is sampled at x0 + i*dx. Requires odd length (even number of intervals).
pub fn simpson_uniform(y: &[f64], dx: f64) -> Result<f64, String> {
    let n = y.len();
    if n < 3 || n % 2 == 0 {
        return Err("simpson_uniform: need odd number of points >= 3".into());
    }
    let mut sum = y[0] + y[n - 1];
    for i in 1..(n - 1) {
        sum += if i % 2 == 0 { 2.0 * y[i] } else { 4.0 * y[i] };
    }
    Ok(sum * dx / 3.0)
}
