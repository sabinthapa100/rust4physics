# rust4physics

**Sabin Thapa** (Kent State University, PhD student) — `sthapa3@kent.edu`  
A small, fast, reproducible **physics lab in Rust**, motivated by quarkonium / Schrödinger numerics (and inspired by QTraj-style workflows).

This repository is structured like a “senior-scientist codebase”:

- `src/lib.rs` = reusable **toolbox** (numerics + quantum utilities)
- `src/bin/*.rs` = reproducible **experiments** (each makes plots/data)
- `cargo test` = physics sanity checks (Hermiticity, trace, etc.)
- outputs go to `out/` by default

---

## Why your HO wavefunction tails oscillate (and how we fix it)

Your attached plot shows oscillatory tails far outside the classical turning point.  
That’s a classic **shooting-method numerical issue**: integrating from one boundary picks up a tiny admixture of the *growing* solution which eventually dominates.

**Fix used here:** integrate from the **origin** with correct **parity** (even/odd), apply **renormalization** during integration to avoid overflow, and use a boundary residual to bracket energies.

You’ll get a clean Gaussian-like ground state.

---

## Quickstart

### 0) Check Rust toolchain
```bash
rustc --version
cargo --version
```

### 1) Run tests (sanity checks)
```bash
cargo test
```

### 2) Harmonic oscillator eigenvalues + wavefunction (Numerov)
Writes:
- `out/ho/energies.txt`
- `out/ho/ground.csv`
- `out/ho/ground.png` (ψ and V)

```bash
cargo run --release --bin tise_ho_numerov -- --outdir out/ho
```

### 3) Free Gaussian wavepacket frames (analytic) + video
Writes PNG frames to `out/wp_frames/`.

```bash
cargo run --release --bin wavepacket_free -- --outdir out/wp_frames --frames 120
```

Convert frames → video (recommended, good compression):
```bash
# MP4 (H.264). CRF 18–28 controls quality/size (lower = higher quality).
ffmpeg -framerate 30 -i out/wp_frames/frame_%04d.png \
  -c:v libx264 -pix_fmt yuv420p -crf 22 -preset medium out/wavepacket.mp4
```

Convert frames → GIF (small demos; larger files than mp4):
```bash
magick -delay 3 -loop 0 out/wp_frames/frame_*.png out/wavepacket.gif
```

---

## Quarkonium / Cornell potential (vacuum + screened medium)

We include a **radial Schrödinger solver** (Numerov + shooting) for central potentials.

Potential options include:
- Cornell: `V(r) = -α/r + σ r + C`
- Screened Cornell (toy Debye screening): uses temperature-dependent screening scale μ_D(T)

### Run: Cornell bound states (S, P, …) for bottomonium-like parameters
Writes:
- `out/cornell/levels.txt`
- `out/cornell/state_l0_n1.png`, etc.

```bash
cargo run --release --bin quarkonium_cornell -- --outdir out/cornell
```

### Run: “Bjorken medium scan” (toy)
Evolves T(τ) with Bjorken scaling `T(τ) = T0 (τ0/τ)^{1/3}` and reports binding vs τ.  
This is not a full open quantum system; it’s a **fast, inspectable stepping stone**.

```bash
cargo run --release --bin bjorken_screened_scan -- --outdir out/bjorken
```

---

## What Rust gives you vs C++ and Python (for this kind of work)

### Rust vs C++
- **Memory safety without GC**: eliminates whole classes of bugs (use-after-free, data races).
- **Fearless concurrency**: parallel loops are far safer to introduce.
- **Cargo**: reproducible dependency management + builds + tests + benchmarks.
- Similar performance to C++ when written idiomatically.

### Rust vs Python
- 10×–100× speedups for tight loops **without** needing C extensions.
- Strong types + ownership = fewer silent bugs in numerical code.
- But: Python is still unbeatable for exploratory notebooks; Rust is best when you want a *robust instrument*.

---

## Physics sanity checks included

Run:
```bash
cargo test
```

Tests include:
- **Hermiticity** of a finite-difference Hamiltonian
- **Trace preservation** under unitary similarity transform `ρ -> U ρ U†`

These are “guardrails” so the code stays correct as you extend it.

---

## Next upgrades (great “fun + research-adjacent” roadmap)

1) Add **split-operator time evolution** for wavepacket scattering on Cornell/screened potentials  
2) Add **complex potentials** (Im V) and monitor norm loss, then restore trace using Lindblad terms  
3) Add **density matrix evolution** in a truncated basis (few bound states) with time-dependent medium rates  
4) Add GPU kernels (later) using `wgpu` or by offloading heavy linear algebra

---

## The cargo bin error you saw

`error: no bin target named tise_ho_numerov` means Cargo cannot *see* the file.

Cargo discovers bins from:
- `src/main.rs` (default bin named `rust4physics`)
- `src/bin/*.rs` (each file = one bin target)

So ensure:
```bash
ls -l src/bin
# should include: tise_ho_numerov.rs, wavepacket_free.rs, ...
```

If the file is missing or named differently, Cargo won’t list it.

---

## Author

Sabin Thapa  
Kent State University (Physics PhD)  
Email: `sthapa3@kent.edu`
# rust4physics
