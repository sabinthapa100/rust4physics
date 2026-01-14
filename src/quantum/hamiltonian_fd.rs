/// Build a 1D finite-difference Hamiltonian H = -1/(2m) d2/dx2 + V(x)
/// with Dirichlet boundary conditions on a uniform grid.
pub fn hamiltonian_fd_1d<V>(x0: f64, dx: f64, n: usize, mass: f64, v: V) -> Vec<Vec<f64>>
where
    V: Fn(f64) -> f64,
{
    let mut h = vec![vec![0.0_f64; n]; n];
    let kin = 1.0 / (2.0 * mass * dx * dx);

    for i in 0..n {
        let x = x0 + i as f64 * dx;
        h[i][i] = 2.0 * kin + v(x);
        if i + 1 < n {
            h[i][i + 1] = -kin;
            h[i + 1][i] = -kin;
        }
    }
    h
}

/// Check symmetry (Hermiticity for real matrices) within tolerance.
pub fn is_symmetric(h: &[Vec<f64>], tol: f64) -> bool {
    let n = h.len();
    for i in 0..n {
        for j in (i + 1)..n {
            if (h[i][j] - h[j][i]).abs() > tol {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quantum::potentials::harmonic_oscillator;

    #[test]
    fn fd_hamiltonian_is_symmetric() {
        let n = 200usize;
        let x0 = -5.0;
        let x1 = 5.0;
        let dx = (x1 - x0) / (n as f64 - 1.0);

        let h = hamiltonian_fd_1d(x0, dx, n, 1.0, |x| harmonic_oscillator(x, 1.0));
        assert!(is_symmetric(&h, 1e-12));
    }
}
