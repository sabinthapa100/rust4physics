/// Simple Bjorken cooling law (ideal, 1D boost-invariant):
///   T(τ) = T0 * (τ0/τ)^(1/3)
pub fn bjorken_temperature(tau: f64, tau0: f64, t0: f64) -> f64 {
    t0 * (tau0 / tau).powf(1.0 / 3.0)
}

/// Toy Debye screening scale mu(T) = a*T
/// (In a more faithful model: m_D ~ g(T)*T * sqrt(1 + Nf/6), etc.)
pub fn mu_debye(t: f64, a: f64) -> f64 {
    a * t
}
