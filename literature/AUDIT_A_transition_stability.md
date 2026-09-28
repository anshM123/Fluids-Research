# Literature death audit A — transition & stability (returned 2026-09-28)

(Agent audit; references marked * were not re-verified live. WebFetch to publishers/arXiv blocked; WebSearch works.)

| Idea | Verdict | Closest work | What remains open |
|---|---|---|---|
| IDEA001 2D PPF transient chaos / 2D puffs / DP | **KNOWN (core)** | Jiménez 1990 JFM 218 (permanent localized patches); Xiao, Tao & Zhang 2021 PoF 33:031706 (localized wave packets, splitting above Re≈4950, "first-order"); Zhang & Tao arXiv:2204.01913 (lifetime ∝ (Re_c−Re)^{-1/2}, deterministic fold); Huang et al. 2024 JFM 994:A6 (2D turbulence sustained Re≳2400); Camobreco, Pothérat & Sheard 2023 JFM 963:R2 (quasi-2D subcritical) | Noise-driven escape vs fold; quasi-2D (Hartmann) DP test |
| IDEA003 sharp 2D Couette threshold γ | **OPEN (numerically)** | Bedrossian–Vicol–Wang 2018 (≤1/2); Masmoudi–Zhao 2022 (≤1/3); Chen–Li–Wei–Zhang 2020 (channel ≤1/2); Vanneste–Morrison–Warn 1998 (echo transient growth) | No numerical estimate of sharp Sobolev γ |
| IDEA057 Kolmogorov+friction transition order | PARTIAL | Bouchet & Simonnet 2009; Frishman–Laurie–Falkovich 2017; van Kan & Alexakis 2019; Novotný et al. 2025 PRFluids 10:054605 (arXiv:2406.08566, experiment, superfluid He, multistability & critical behaviour) | Strict 2D (Re, α, aspect) transition-order map |
| IDEA021 deflation in canonical benchmarks | PARTIAL | Farrell–Birkisson–Funke 2015; Boullé–Dallas–Farrell 2022 (RBC); Auteri et al. 2002 (cavity Hopf ≈ 8018); multi-lid multiplicity known | No deflation search of 1-lid cavity / single cylinder; stable new branches unlikely |
| IDEA053 minimal-seed universality | PARTIAL (universality likely false) | Pringle–Kerswell 2010; Duguet et al. 2013; Zhang et al. 2023 PoF 35:051704 (2D PPF E_c ∝ Re^{-3.8}) | — |
| IDEA052 super-exponential lifetimes mechanism | PARTIAL, contested | Goldenfeld–Guttenberg–Gioia 2010; Hof 2006/08; Linkmann–Morozov 2015 PRL; Kreilos–Eckhardt–Schneider 2014; Shih–Hsieh–Goldenfeld 2016 Nat Phys; Nemoto–Alexakis 2021; Gomé–Tuckerman–Barkley 2022; Guan & Tao 2025 arXiv:2504.14465 | Controlled cross-PDE discrimination of mechanisms (rare-event sampling) |
| IDEA063 2D cylinder drag crisis | PARTIAL | Singh & Mittal 2005 (2D FEM crisis at 2×10⁵); Tamura et al. 1990 | grid-converged, dissipation-free 2D at 10⁶ — beyond budget |
| IDEA010 tripole lifetime | KNOWN small amp / PARTIAL large | Bernoff & Lingevitch 1994; Bajer–Bassom–Gilbert 2001; Rossi–Lingevitch–Bernoff 1997; Gallay 2018 | large-amplitude threshold ε_c(Re) |
| IDEA005 merger critical slowing | **OPEN (as a law)** | Dritschel 1985/1995; Melander–Zabusky–McWilliams 1988; Meunier et al. 2002 (a/b_c≈0.24); Josserand & Rossi 2007 | No measured (d_c−d)^{-1/2}; slow-passage prediction (a/b)_c(Re)−(a/b)_c(∞) ∝ Re^{-2/3}, lag ∝ Re^{1/3} untested |
| IDEA041 non-shear steady states near Kolmogorov/Poiseuille | OPEN (numerically) | Lin & Zeng 2011; Coti Zelati–Elgindi–Widmayer 2023; Castro & Lear 2023; arXiv:2404.19034; arXiv:2609.18280 | no numerical continuation of these families |

Additional open questions surfaced: noise on 2D PPF wave packet (fold → Kramers); Couette-fraction threshold for sustained 2D chaos (Couette–Poiseuille mix; Falkovich & Vladimirova 2018 PRL: pure 2D Couette cannot sustain); switching between symmetric/asymmetric 2D channel turbulence (Markeviciute & Kerswell 2021 JFM 917:A57); wave packet as localized RPO & snaking; quasi-2D duct DP test.

**Decisions:** IDEA001 KILLED (core known) → second-generation variants IDEA101 (Couette–Poiseuille fraction threshold for sustained 2D chaos) and IDEA102 (noise-activated escape of 2D wave packets) logged. IDEA005 & IDEA003 promoted to cheap tests. IDEA053 KILLED (universality likely false; negative-result trajectory). IDEA063 KILLED (budget). IDEA021 → Tier C.
