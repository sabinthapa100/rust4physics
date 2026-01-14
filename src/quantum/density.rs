use num_complex::Complex64;

/// Construct pure-state density matrix ρ = |ψ><ψ| from a state vector ψ.
pub fn density_from_state(psi: &[Complex64]) -> Vec<Vec<Complex64>> {
    let n = psi.len();
    let mut rho = vec![vec![Complex64::new(0.0, 0.0); n]; n];
    for i in 0..n {
        for j in 0..n {
            rho[i][j] = psi[i] * psi[j].conj();
        }
    }
    rho
}

/// Trace of a density matrix.
pub fn trace(rho: &[Vec<Complex64>]) -> Complex64 {
    let n = rho.len();
    let mut tr = Complex64::new(0.0, 0.0);
    for i in 0..n {
        tr += rho[i][i];
    }
    tr
}

/// Apply similarity transform: ρ -> U ρ U†
pub fn unitary_similarity(u: &[Vec<Complex64>], rho: &[Vec<Complex64>]) -> Vec<Vec<Complex64>> {
    let n = rho.len();
    let mut tmp = vec![vec![Complex64::new(0.0, 0.0); n]; n];
    let mut out = vec![vec![Complex64::new(0.0, 0.0); n]; n];

    // tmp = U * rho
    for i in 0..n {
        for k in 0..n {
            let uik = u[i][k];
            if uik == Complex64::new(0.0, 0.0) { continue; }
            for j in 0..n {
                tmp[i][j] += uik * rho[k][j];
            }
        }
    }

    // out = tmp * U†
    for i in 0..n {
        for j in 0..n {
            let mut s = Complex64::new(0.0, 0.0);
            for k in 0..n {
                s += tmp[i][k] * u[j][k].conj(); // (U†)_{k j} = conj(U_{j k})
            }
            out[i][j] = s;
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn trace_preserved_under_unitary_similarity() {
        // 2D pure state
        let psi = vec![Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0), Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)];
        let rho = density_from_state(&psi);

        // U = rotation in SU(2): [[cosθ, sinθ],[-sinθ, cosθ]]
        let theta: f64 = 0.37;
        let c = theta.cos();
        let s = theta.sin();
        let u = vec![
            vec![Complex64::new(c, 0.0), Complex64::new(s, 0.0)],
            vec![Complex64::new(-s, 0.0), Complex64::new(c, 0.0)],
        ];

        let rho2 = unitary_similarity(&u, &rho);
        let tr1 = trace(&rho);
        let tr2 = trace(&rho2);

        assert_relative_eq!(tr1.re, 1.0, epsilon = 1e-12);
        assert_relative_eq!(tr2.re, 1.0, epsilon = 1e-12);
        assert_relative_eq!(tr2.im, 0.0, epsilon = 1e-12);
    }
}
