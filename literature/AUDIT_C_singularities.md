# Literature death audit C — singularities & mathematical fluid dynamics (returned 2026-09-28)

| Idea | Verdict | Closest | Open / angle |
|---|---|---|---|
| IDEA028 unstable-singularity ladders | PARTIAL; **core question explicitly open** | Wang et al. arXiv:2509.14185 (2025): CCF stable + 1st unstable (λ₁≈0.6057, first found by Eggers–Fontelos, Nonlinearity 33:325, 2020) + NEW 2nd unstable λ₂=0.4703 (residual ~1e-7); empirical laws only for IPM λₙ≈1/(1.1459n+0.9723) and Boussinesq λₙ≈1+1/(1.4187n+1.0863); for CCF "no clear asymptotic relation … premature". Link: dissipative CCF (−Δ)^{α/2} blow-up expected for α ≤ 1/(1+λₙ) → 0.623, 0.68. Wang–Léger–Lai–Buckmaster arXiv:2511.22819: 4th unstable IPM; **CCF 3rd profile not found** (scan λ∈[0.455,0.4713] leaves non-smooth signal at origin). Huang–Qin–Wang–Wei ARMA 2024 (gCLM a≤1 smooth profiles); Huang–Tong–Wei CMP 2023 (De Gregorio countably many, same scaling); Huang–Tong–Wang arXiv:2603.25104 (degenerate data → gCLM unstable ladder); Lushnikov–Silantyev–Siegel JNS 2021 (a_c≈0.6891); Chen–Hou–Huang CPAM 2021; Elgindi–Jeong ARMA 2020; Zheng Nonlinearity 2023; Xu 2026 (arXiv:2607.19762, 2609.13220); Hou–Luo: Chen–Hou–Huang Ann PDE 2022, CMP 2025 (arXiv:2308.01528), Chen–Huang–Li arXiv:2604.01868; Hou group arXiv:2506.19243 | existence of CCF n≥3 smooth profiles; large-n law (λₙ→0 ⇒ α_c→1; accumulation λ*>0 caps it); Hou–Luo 1D unstable ladder |
| IDEA037 CAP chaos in 2D NS | OPEN (moderate confidence) | Wilczak–Zgliczyński JDE 2020 (KS symbolic dynamics); van den Berg–Breden–Lessard–van Veen JNS 2021; Arioli–Koch JMFM 2021 (Hopf planar NS); Bedrossian–Punshon-Smith CMP 2024 (stochastic Galerkin NS N≥392); Wilczak–Zgliczyński JDE 2026 (arXiv:2502.09760) | heavy; first milestone = CAP of a Kolmogorov periodic orbit |
| IDEA039 ∇ω growth on T² | OPEN | Kiselev–Šverák 2014; Zlatoš 2015; Denisov 2009/2015; Jeong–Yao et al. arXiv:2507.15739 | numerics only guide |
| IDEA040 unforced Leray–Hopf non-uniqueness | **KNOWN (R³)** | Hou–Wang–Yang arXiv:2509.25116 (CAP, 2025); 2D numerics Albritton–Guillod–Korobkov–Ren arXiv:2601.03161 (three self-similar solutions σ∈[39.2,80]) | CAP in 2D |
| IDEA042 N-vortex collapse counting | PARTIAL | Novikov–Sedov 1979; Kudela JNS 2014; O'Neil 1987; arXiv:2410.14973; 2505.19782; 2607.16490 | counting law via certified homotopy continuation |
| IDEA044 Boussinesq fractional dissipation | regularity partial; blow-up open | Math Ann 391 (2025); arXiv:2606.03680; Chen–Hou | threshold ladder α_c(n) from λₙ |
| IDEA055 viscously arrested singularity | PARTIAL (classical matching) | Eggers–Fontelos book 2015; Schochet 1986 (viscous CLM blows up) | modest novelty |
| IDEA051 max amplification exponents | PARTIAL | Ayala–Protas 2011/2014; Ayala–Doering–Simon 2018; Yun–Protas 2018; Kang–Yun–Protas 2020 | derivation of α from scaling |
| IDEA043 sharp constants | PARTIAL | Lu–Doering 2008 | certified extremals |

Additional: CAP of 2D NS non-uniqueness (AGKR profiles); Hou–Luo 1D smooth unstable ladder vs Boussinesq law; certified unstable-eigenvalue counts for gCLM vs a; De Gregorio near −sin kθ (k≥3); gSQG point-vortex collapse/burst.
Context (unverified): 8 Sep 2026 manuscripts reportedly claim unforced Euler blow-up from smooth data on R³ and forced NS blow-up.

**Decisions:** IDEA028 → **TIER A (P02)**. IDEA040 KILLED (known). IDEA037 PAUSED (feasibility). IDEA044 MERGED into P02 (threshold ladder). IDEA055 MERGED (modest). IDEA042 Tier C.
