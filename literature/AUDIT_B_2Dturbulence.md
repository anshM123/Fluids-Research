# Literature death audit B — 2D turbulence & long-time Euler (returned 2026-09-28)

| Idea | Verdict | Closest work | Open |
|---|---|---|---|
| IDEA002 integrability selection | PARTIAL | Modin & Viviani JFM 884 A22 (2020) sphere: 4 blobs at L≈0, else 3/2; survey Arnold Math J (2021) states torus version; Modin & Viviani ARMA 2025 (arXiv:2405.14282), J Comput Dyn 2026 (arXiv:2508.07088): blob number correlates with abs(L)/sqrt(enstrophy); Dritschel–Qi–Marston JFM 783 (2015); rotating sphere Ryono et al. arXiv:2606.27778; counter-evidence: vortex crystals in disks (Fine et al. PRL 1995; Jin & Dubin PRL 1998) | no ensemble test on flat torus/disk; no quantitative law for surviving-vortex number off the sphere |
| IDEA009 condensate switching | PARTIAL | Bouchet & Simonnet 2009; Bouchet–Laurie–Zaboronski 2011/2014; Laurie & Bouchet NJP 2015; Xu, van Kan, Liu & Knobloch PRF 9 064605 (2024) Kramers-like DNS lifetimes; Bouchet–Rolland–Simonnet PRL 2019 (β-plane jets, AMS) | computed non-equilibrium quasi-potential reproducing measured rates |
| IDEA011 D_KY(Re) 2D | LARGELY KNOWN | Grappin & Léorat JFM 1991; Platt–Sirovich–Fitzmaurice 1991; "Dimensional regimes in Kolmogorov flow" arXiv:2602.08960 (PRF 2026): D_KY saturates, saturation ∝ n; CFT bound; Liu CMP 1993 | saturation at higher Re? |
| IDEA034 manifold dimension | KNOWN | Page–Brenner–Kerswell PRF 2021; Pérez De Jesús & Graham PRF 2023; PRE 2025 (doi 10.1103/cdz2-858n) | — |
| IDEA016 turbulence on hyperbolic surfaces | **OPEN (numerics)** | Falkovich & Gawędzki JSP 156 (2014) theory; Khesin & Misiołek PNAS 2012; Reuther & Voigt 2015/2018 (curved surfaces, no cascade diagnostics); arXiv:2604.25682 point vortices | first forced-turbulence DNS on H² testing arrest at curvature radius |
| IDEA048 ergodic optimisation | target killed | Hunt & Ott PRL 1996; Yang–Hunt–Ott 2000; Contreras 2016; Tobasco–Goluskin–Doering 2018; Lakshmi et al. SIADS 2020 | for Kolmogorov forcing ⟨D⟩ ≤ D_lam (Cauchy–Schwarz) → max is laminar; non-trivial targets remain |
| IDEA035 UPO statistics | PARTIAL, active | Chandler & Kerswell 2013; Lucas & Kerswell 2015; Page et al. PNAS 2024; Farazmand & Sapsis 2017 | tails/large deviations from UPO weights |
| IDEA004 wall anomalous dissipation | NOT SETTLED | Clercx & van Heijst PRE 2002 (→0); Nguyen van yen et al. PRL 2011 (const, penalization); JFM 2018; Kato 1984; Constantin & Vicol 2018; Bardos–Titi–Wiedemann 2019 | converged boundary-fitted D(Re) at Re ≳ 1e5 |
| IDEA025 Zeitlin torus | PARTIAL | Zeitlin 1991; Abramov & Majda PNAS 2003; Dubinkina & Frank 2007/10; Qi & Marston 2014; Cifani et al. PRF 2022 | merge with IDEA002 |
| IDEA056 truncated vs continuum invariants | PARTIAL | Kraichnan; Fox & Orszag 1973; Basdevant & Sadourny 1975; Abramov & Majda 2003; Bouchet & Corvellec 2010; Ray et al. PRE 2011 | crossover time t*(N) law |

Further open questions surfaced: viscous lifetime of Modin–Viviani blob states; continuation of non-shear Euler states near Kolmogorov (dup. IDEA041); Navier-slip ℓ_s ∝ ν^γ threshold for anomalous wall dissipation; long-time Euler on a torus of revolution (Sakajo & Shimizu 2016); k^{-5} large-scale range in drag-free condensates (van Kan, Alexakis & Knobloch arXiv:2504.02978).

**Decisions:** IDEA011, IDEA034 KILLED (known). IDEA048 KILLED as stated. IDEA016 → Tier B candidate (open, conceptually clean, but geometry/solver cost; compact-surface caveat). IDEA002 → Tier B (torus/disk ensemble law). IDEA004 → Tier C (Chebyshev channel solver exists; high-Re resolution cost). IDEA009 → PAUSED (Bouchet group territory).
