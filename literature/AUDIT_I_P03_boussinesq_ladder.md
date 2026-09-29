# AUDIT I — novelty audit for P03 (2D Boussinesq ladder, WKB mechanism)

Environment limits: arxiv.org and most publisher, news and blog domains are blocked. The audit therefore relies
on search-engine snippets, abstracts, and the repository's earlier audits (AUDIT_C, AUDIT_H).

| Claim | Nearest literature | Verdict |
|---|---|---|
| C1. All known smooth 2D Boussinesq self-similar profiles (stable + unstable) lie on ONE branch of least-singular profiles; they are the zeros of m(λ) − 2. | Wang et al. arXiv:2509.14185 (PINN discovery of λ₀…λ₃ + candidate; each profile found separately; empirical law). Chen–Hou (stable profile; arXiv:2210.07191). Chen–Huang–Li arXiv:2604.01868 (singular, unbounded profiles for Hou–Luo and Boussinesq, a different class). | **New.** No branch or continuation picture for Boussinesq is reported. The least-singular family idea parallels our CCF work (P02), and families of non-smooth profiles exist for gCLM (Huang–Tong–Wang arXiv:2603.25104). Credit both. |
| C2. Classical computation of λ₀…λ₇ (≥ 8 smooth profiles); λ₅–λ₇ new. | Wang et al.: 3 unstable + 1 candidate. The follow-up arXiv:2511.22819 concerns IPM/CCF. | **New** (as far as accessible). Their tabulated λ₂, λ₃ were not accessible. Ours differ by 2.5e-3 and 2.9e-3 from their two-point law, which passes through λ₀ and λ₁. |
| C3. As λ → 1: quasi-stagnant boundary region plus a front tending to a square-root cusp; no sonic point at λ > 1. | CCF terminal cusp (our P02). No Boussinesq analogue found. | **New.** |
| C4. Ladder = WKB (phase) quantisation in the stalled (quasi-stagnant) layer; ΔRe Φ → π; asymptotic spacing π/C between 1.476 and ≈ 1.51 (1.478 with analytic corrections, 1.49–1.50 with the dip-induced non-analytic corrections; 3/2 a candidate); smoothness defect ∝ (λ−1)^{p} e^{−Im Φ} with fitted p ≈ 1.5. Formal asymptotics + an elementary lemma; no proof claimed. | Exponential-asymptotics selection (Saffman–Taylor, dendrites; Kruskal–Segur type) is classical in other contexts. Infinitely many self-similar blow-up profiles are known for Keller–Segel in d = 3–9 (arXiv:2503.02263, matched asymptotics). De Gregorio: infinitely many profiles with the same scaling (Huang–Tong–Wei, CMP 2023). | **New for fluid singularities.** Frame it as a WKB/beyond-all-orders selection mechanism and cite the classical selection literature. |
| C5. Hou–Luo model: the same ladder (12 crossings), λ₀^{HL} = 1.99871. | Huang–Qin–Wang–Wei CMP 2025 (existence of the smooth stable HL profile); Chen–Hou–Huang Ann. PDE 2022. The HL unstable ladder was listed as open in our AUDIT_C. | **New.** Check λ₀^{HL} against the published stable exponent when accessible. |
| C6. Finite (CCF) vs infinite (Boussinesq/HL) ladders explained by where the sonic cusp forms. | — | **New synthesis.** |

**Kill check:** no source found that computes the Boussinesq unstable ladder beyond the four PINN profiles, that
establishes its infiniteness, or that explains the empirical law. Status: **survives**. Before submission, re-check
arXiv once access is available, specifically: Wang et al. v2+, 2511.22819, any 2026 work on "unstable Boussinesq
profiles", and Hou group papers.

**Addendum (R028, reviewer alignment).** Three further claims.

| Claim | Nearest literature | Verdict |
|---|---|---|
| C7. Instability index n of the n-th profile for all eight profiles, from two independent linear-stability methods (march-based eigen-condition; global shift-invert eigen-solver); 21 unstable eigenvalues agree to 1e-5. | Wang et al. 2509.14185 report the n-unstable-modes pattern for their profiles (n ≤ 3). | **Extends** their observation to n = 4–7. The two classical stability methods are new for this problem. |
| C8. Independent reproduction of λ₀…λ₇ by a second solver sharing no numerical ingredient, to ≤ 1e-5. | — | Verification (no novelty claim). |
| C9. The lower unstable eigenvalues scale with ε = (λ_n−1)/2 (μ_k/ε → 1.6, 3.8, 6.05, …). | — | **Observation**, not derived. Check for related eigenvalue asymptotics in the literature before claiming. |
