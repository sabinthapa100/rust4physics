/// Harmonic oscillator potential V(x) = 0.5 * ω^2 x^2 (ħ=1, m=1 conventions elsewhere).
pub fn harmonic_oscillator(x: f64, omega: f64) -> f64 {
    0.5 * omega * omega * x * x
}

/// Square barrier potential: V(x)=v0 on [a,b], else 0.
pub fn square_barrier(x: f64, a: f64, b: f64, v0: f64) -> f64 {
    if x >= a && x <= b {
        v0
    } else {
        0.0
    }
}

/// Cornell potential: V(r) = -alpha/r + sigma*r + c
/// r in [GeV^-1] if using natural units, or [fm] if you consistently convert.
pub fn cornell(r: f64, alpha: f64, sigma: f64, c: f64) -> f64 {
    let eps = 1e-12;
    -alpha / (r.max(eps)) + sigma * r + c
}

/// Toy screened Cornell (Debye-screened) potential.
/// A commonly used smooth form is:
///   V(r,T) = -alpha * exp(-mu r)/r + (sigma/mu) * (1 - exp(-mu r)) + c
/// where mu ~ m_D(T) is screening scale.
///
/// This is a *toy* (fast, inspectable). For production physics you may tune form/params.
pub fn screened_cornell(r: f64, alpha: f64, sigma: f64, mu: f64, c: f64) -> f64 {
    let eps = 1e-12;
    let rr = r.max(eps);
    let e = (-mu * rr).exp();
    -alpha * e / rr + (sigma / mu) * (1.0 - e) + c
}
