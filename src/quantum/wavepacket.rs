use std::f64::consts::PI;

/// Free-particle Gaussian wavepacket (analytic) in 1D, units ħ=1, m=1.
///
/// We return probability density |psi(x,t)|^2:
///   |psi|^2 = (1 / sqrt(2π σ_t^2)) * exp(-(x-x_c(t))^2 / (2 σ_t^2))
/// with σ_t^2 = σ0^2 (1 + (t/(2σ0^2))^2), and x_c(t)=x0 + k0 t
pub fn free_gaussian_density(x: f64, t: f64, x0: f64, sigma0: f64, k0: f64) -> f64 {
    let xc = x0 + k0 * t;
    let tau = t / (2.0 * sigma0 * sigma0);
    let sigma_t2 = sigma0 * sigma0 * (1.0 + tau * tau);

    let norm = 1.0 / ( (2.0 * PI * sigma_t2).sqrt() );
    let arg = -((x - xc) * (x - xc)) / (2.0 * sigma_t2);

    norm * arg.exp()
}
