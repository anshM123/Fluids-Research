# Second literature death audit — P02 main claims (2026-09-28)

Search budget was nearly exhausted after audits A–G; arxiv.org, ar5iv and pith.science are blocked by the
egress proxy, so verdicts rest on search-engine abstracts plus the earlier audit notes (AUDIT_C). Queries
run for this audit: "Córdoba–Córdoba–Fontelos third unstable self-similar profile finite number";
"self-similar blow-up profiles branch spiral log-periodic sonic point square-root singularity nonlocal
transport Hilbert transform"; "Córdoba–Córdoba–Fontelos self-similar blow-up 2026 unstable profiles
continuation classical numerical"; plus targeted look-ups of the hits below.

| Claim | Closest work found | Verdict |
|---|---|---|
| C1. All known smooth CCF profiles (λ₀ stable, λ₁, λ₂) lie on ONE connected branch of self-similar profiles; smooth profiles = crossings p(λ)=2 of the local exponent at the origin; instability index increases by one per crossing | Eggers–Fontelos, Nonlinearity 33 (2020) (λ₁); Wang et al. arXiv:2509.14185 (λ₂=0.4703, empirical λ_n laws for IPM/Boussinesq, "premature" for CCF); Wang–Léger–Lai–Buckmaster arXiv:2511.22819 (λ₂=0.47132422; λ₃ scan in [0.455, 0.4713] finds only a non-smooth signal at the origin). Continuous families of *singular* profiles (non-integer vanishing order at the origin) are used for gCLM with degenerate data: Huang–Tong–Wang arXiv:2603.25104 (2026); Hou–Luo/Boussinesq singular profiles: Chen–Huang–Li arXiv:2604.01868 (2026) | **OPEN / new** for CCF (global connectivity + crossing mechanism not reported). The notion of a λ-family with non-smooth origin behaviour is known in related models → must be credited. |
| C2. The branch terminates (sonic depth δ→0) at a *cusp profile* with an interior square-root singularity, B² = 2λΘ_s, layer width w ∝ δ² | Hoang–Radosz, ARMA 2017 (arXiv:1602.02451): time-dependent cusp/needle formation for a CCF-*inspired* nonlocal active scalar | **OPEN / new** (different object: a terminal self-similar profile, not a time-dependent cusp; cite as related). |
| C3. Log-periodic (spiral) approach to the terminal cusp with universal frequency 2τ, τ tanh(πτ/2) = 1/2, from the local operator u' = −H[u]/(2|x|) | Classical spirals of solution branches approaching singular solutions (Joseph–Lundgren 1973 Gelfand problem; Emden–Fowler); no Hilbert-transform/nonlocal analogue found | **OPEN / new** in this context (mechanism analogous to classical spirals; cite). |
| C4. The CCF principal branch carries exactly three smooth self-similar profiles (no λ₃ on it); p stays in (2.0050, 2.0134) on the whole arc beyond λ₂ | arXiv:2511.22819 search failure in [0.455, 0.4713] (consistent) | **OPEN / new** (answers the question left open by 2509.14185/2511.22819 for the principal branch; does not exclude disconnected branches). |
| C5. Classical log-variable spectral continuation reproduces λ₂ to 12 digits in ≈20 s on one CPU core, where PINN+Gauss–Newton pipelines needed specialised gradient-normalised losses | arXiv:2511.22819 (PINN machine precision), 2509.14185 | Methodological comparison; fair only for 1D. |

Additional related work to cite: Elgindi–Ghoul–Masmoudi, Anal. PDE 14 (2021) (stable self-similar blow-up for nonlocal transport/stretching families); Chen–Hou–Huang (Hou–Luo, gCLM); Lushnikov–Silantyev–Siegel JNS 2021 (gCLM a_c≈0.6891); Silvestre–Vicol 2016 (CCF/transport with nonlocal velocity); Córdoba–Córdoba–Fontelos, Ann. Math. 2005.

**Decision:** no kill. Claims C1–C4 proceed to the dossier with the scope restriction "on the connected branch containing all known profiles".
