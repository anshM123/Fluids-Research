# Literature death audit E — computational methods & scientific ML (returned 2026-09-28)

| Idea | Verdict | Closest | Open / angle |
|---|---|---|---|
| IDEA022 tensor-network complexity | PARTIAL | Gourianov et al. NCS 2022 (χ saturates in 2D); Hölscher et al. PRR 7 013112 (2025): χ=O(poly(1/ε)) from spectrum, plateau for Re up to 1e7 (2D); Pisoni et al. PRR (3D TT-DNS Re_λ=315); Horner et al. arXiv:2606.17064 (1D exact when χ ≥ #active modes); 2D RB arXiv:2604.16179 (snapshot χ grows to Ra=1e11); Esmaeili…Jaksch arXiv:2608.28869 | volume law χ vs domain size / number of eddies; KS & forced Burgers tests; statistics-χ vs snapshot-χ; honest cost crossover N ≫ χ³–χ⁴ |
| IDEA023 time-parallel statistics | PARTIAL (≈known ODE) | Wang et al. PoF 2013 (LSS); Samaddar/Reynolds-Barredo JCP 2010–13; Lunet et al. CVS 2018; Vargas et al. SISC 2023; arXiv:2604.00855 | trivial ensemble baseline strong → KILL |
| IDEA026 topology-certified AMR | OPEN in solvers | cpSZ/TspSZ compression; Finken et al. arXiv:2608.12142; Chen–Mischaikow–Laramee–Zhang TVCG 2008; Kamkar et al. JCP 2011 | per-cell topological-degree certificates; Moffatt eddy benchmark |
| IDEA027 certified shadowing for PDE turbulence | OPEN for PDEs | Hammel–Yorke–Grebogi; Coomes–Koçak–Palmer 1995; Hayes & Jackson 2003/05; Larsson & Sanz-Serna 1999; Chandramoorthy & Wang 2021; Budanur 2024 | computable trust horizon for KS/2D NS incl. truncation; could arbitrate Qin & Liao (CNS, JFM 2022) vs McMullen et al. arXiv:2510.04828 |
| IDEA029 walk-on-spheres Stokes | little in fluids proper | Rioux-Lavoie et al. TOG 2022; Sugimoto–Batty–Hachisuka 2024; Walk on Stars (Sawhney 2023); NIST ZENO | unbiased grid-free no-slip Stokes estimator (Papkovich–Neuber / walk-on-boundary) |
| IDEA030 exp. integrators for viscoelastic | premise flawed | Gupta & Vincenzi 2019; Beneitez–Page–Kerswell 2023; Yerasi et al. 2024 | pivot: flow-map (characteristic-mapping) Lodge-integral solver, positive-definite by construction (check Hulsen deformation-fields method) |
| IDEA031 sync threshold law | PARTIAL | Yoshida et al. 2005; Lalescu et al. 2013; Clark Di Leoni et al. PRX 2020; Inubushi et al. PRL 2023; Olson–Titi; Azouani–Olson–Titi 2014 | first-principles k_cη criterion across KS/2D/Sabra |
| IDEA032 hidden-symmetry closure | PARTIAL | Biferale–Mailybaev–Parisi PRE 2017; Mailybaev 2021/22; Campolina & Mailybaev log-lattices | NS LES closure with anomalous exponents |
| IDEA047 QD regime discovery | OPEN for regimes | SAIL 2018; Grizou et al. Sci Adv 2020; Etcheverry et al. 2020 | MAP-Elites over forcing space with invariant descriptors |
| IDEA049 dimensionless learning | LARGELY KNOWN | Bakarji et al. NCS 2022 (BuckiNet); Xie et al. Nat Commun 2022 | incomplete similarity / anomalous exponents |

Additional directions: MLMC couplings for chaotic statistics; shadowing-certified multiple shooting; hidden-symmetry renormalized solvers on log-lattices.

**Decisions:** IDEA023 KILLED; IDEA049 KILLED; IDEA030 PIVOT → IDEA103 (flow-map Lodge-integral viscoelastic solver) Tier C; IDEA022 → Tier B pivot: "volume law" test (IDEA104) — but outcome likely a limitation (negative for TN) → keep cheap only; IDEA027 → Tier B candidate (certified trust horizon; could arbitrate CNS dispute); IDEA029 → Tier B candidate (grid-free no-slip Stokes; high method upside, high risk).
