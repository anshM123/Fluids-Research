# Literature death audit D — mixing, transport, low-Re (returned 2026-09-28)

| Idea | Verdict | Closest | Open / angle |
|---|---|---|---|
| IDEA038 mixing speed limit | PARTIAL, fast-moving | Lin–Thiffeault–Doering 2011; Seis 2013; Iyer–Kiselev–Xu 2014; Alberti–Crippa–Mazzucato 2019; Elgindi–Liss–Mattingly 2023; arXiv:2608.14346 (shell model global opt.); Hu–Li–Zhang–Zuazua arXiv:2601.06294; Li–Zuazua arXiv:2605.04688 | sharp c* in ‖θ‖_{H⁻¹} ≳ e^{−c*Γt/L}; Floquet formulation |
| IDEA012 Taylor-dispersion isoperimetry | **FALSE as posed** (thin rectangles beat disk at fixed area & mean velocity; K/(U²A)≈0.038/AR→0) | Doshi–Daiya–Gill 1978; Chatwin–Sullivan 1982; Ajdari–Bontoux–Stone 2006; Aminian et al. 2015/2016 | agent FEM: disk strict local min at fixed flux; **disk global MAXIMISER at fixed pressure gradient** (n-gons non-monotone) — Talenti-type inequality candidate |
| IDEA058 N Stokeslets escape rate | PARTIAL; N-scaling OPEN | Jánosi et al. 1997; Ekiel-Jeżewska et al. 2008; Gruca et al. 2015; arXiv:2501.04965 | κ(N), super-transients τ~e^{αN} |
| IDEA064 roughness shape optimisation | PARTIAL | Toppaladoddi et al. 2015/2017; Zhu et al.; Goluskin–Doering 2016; Wen–Goluskin–Doering 2022 | adjoint free-form wall for steady rolls vs Ra^{1/2} |
| IDEA066 max Lagrangian chaos per enstrophy | PARTIAL | D'Alessandro–Dahleh–Mezić 1999; Finn–Thiffeault 2011; Khakhar–Rising–Ottino 1986; fastest dynamo (Alexakis 2011; Willis 2012) | R*∈[0.36, 0.5] (agent computations: sine shears 0.27, tent shears 0.36; bound 1/2) |
| IDEA013 2D caustics | PARTIAL | Wilkinson–Mehlig; Gustavsson–Mehlig 2016; Meibohm et al. 2021 | C in 2D |
| IDEA059 finger-width principle | KNOWN → KILL | McLean–Saffman; Combescot et al.; Tanveer; Chapman 1999 | — |
| IDEA071 wall-to-wall transport | MOSTLY KNOWN | Tobasco–Doering 2017/2019; Souza–Tobasco–Doering 2020; arXiv:2507.12564 | log gap only |
| IDEA072 enhanced dissipation | KNOWN/PARTIAL | Bedrossian–Coti Zelati 2017; Albritton–Beekie–Novack 2022 | — |

**Decisions:** IDEA012 → pivot IDEA108 (fixed-pressure-gradient disk-maximiser inequality) Tier C (math, modest). IDEA066 → Tier B-/C (clean extremal constant, cheap maps). IDEA038 PAUSED (fast-moving, crowded 2026). IDEA058 Tier C. IDEA059/071/072 KILLED.
