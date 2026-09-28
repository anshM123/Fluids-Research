import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from ccf_nk import CCFNK
files = [('ladder_F16k_up_cross2.npy',16384,'stable λ0'),('ladder_F16k_up_cross1.npy',16384,'λ1'),
         ('ladder_F64k_dn_cross1.npy',65536,'λ2'),('ladder_F16k_dn_cross2.npy',16384,'λ3 cand. (under-resolved)')]
fig, ax = plt.subplots(1,3, figsize=(15,4.2))
for f,N,lab in files:
    S = CCFNK(30,120,N,0.7); d=np.load(f); lam=d[0]; phi=d[1:]
    F,den = S.F(phi,lam,S.normval(phi,lam)); xi=np.exp(S.eta); Th=np.exp(phi+S.c*S.eta)
    beta=lam/(1+lam); C=Th[S.iR]/xi[S.iR]**beta
    m=(S.eta>-6)&(S.eta<6)
    ax[0].plot(xi[m], Th[m]/C, label=f'{lab}: λ={lam:.5f}')
    ax[1].semilogx(xi[m], den[m], label=lab)
    ax[2].semilogx(xi[m], np.gradient(np.log(Th[m]),S.eta[m]), label=lab)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('ξ'); ax[0].set_ylabel('Θ/C'); ax[0].legend(fontsize=7)
ax[1].set_xlabel('ξ'); ax[1].set_ylabel('1+λ+HΘ/ξ  (sonic factor)'); ax[1].set_ylim(0,2)
ax[2].set_xlabel('ξ'); ax[2].set_ylabel('d lnΘ / d lnξ'); ax[2].set_ylim(0,12)
plt.tight_layout(); plt.savefig('fig_profiles.png', dpi=120)
print('saved')
