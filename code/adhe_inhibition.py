import sys; sys.path.insert(0,'/home/claude')
import numpy as np, pandas as pd, json
from lib2 import (load_model, make_configured, P, stress_S, kd_func, atpm_burden,
                  r_func, pfl_inhibition, set_coupling)
def S_onoo(t,alpha): return P.w_ONOO*alpha*P.ONOO_0*np.exp(-P.k_ONOO*t)
def I_adhe(t,alpha,K,n=2.0):
    s=S_onoo(t,alpha); return 1.0/(1.0+(s/K)**n)
VMAX0=17.0
def run_adhe(m,alpha,K,n=2.0):
    dt=P.dt;nstep=int(P.T_batch/dt)
    pfl0=m.reactions.PFL.upper_bound
    atpm0=getattr(m,'_atpm0',None)
    if atpm0 is None: atpm0=m.reactions.ATPM.lower_bound;m._atpm0=atpm0
    e0=m.reactions.EX_etoh_e.upper_bound
    X=P.X_0;glc=P.glc_0;etoh=0.0
    for i in range(nstep):
        t=i*dt;S=stress_S(t,alpha)
        set_coupling(m,r_func(S))
        m.reactions.PFL.upper_bound=pfl0*pfl_inhibition(S)
        m.reactions.ATPM.lower_bound=atpm0+atpm_burden(S)
        m.reactions.EX_etoh_e.upper_bound=VMAX0*I_adhe(t,alpha,K,n)
        m.reactions.EX_glc__D_e.lower_bound=0.0 if glc<0.1 else -P.glc_uptake_max
        try:
            sol=m.optimize()
            if sol.status!='optimal': X=X*np.exp(-kd_func(S)*dt); continue
            X=X*np.exp((sol.objective_value-kd_func(S))*dt)
            glc=max(0,glc+sol.fluxes['EX_glc__D_e']*X*dt)
            etoh+=sol.fluxes['EX_etoh_e']*X*dt
        except Exception: pass
    m.reactions.PFL.upper_bound=pfl0;m.reactions.ATPM.lower_bound=atpm0;m.reactions.EX_etoh_e.upper_bound=e0
    return etoh

GSMM=sys.argv[1]; path='/mnt/user-data/uploads/%s.xml'%GSMM
m=load_model(path); cfg=make_configured(m,'WT')
alphas=np.round(np.arange(0,1.51,0.1),2)
Ks=[('none',1e9),('weak',3.0),('moderate',1.5),('strong',1.0),('severe',0.5)]
rows=[]
for label,K in Ks:
    et=np.array([run_adhe(cfg,a,K) for a in alphas])
    base=et[0]; fold=et.max()/base if base>0 else 0; astar=float(alphas[int(np.argmax(et))])
    rows.append({'inhibition':label,'K_AdhE':K,'GSMM':GSMM,'baseline':round(base,2),
                 'peak':round(et.max(),2),'fold':round(fold,3),'alpha_star':astar,
                 'I_op':round(I_adhe(0,0.6,K),3),'curve':[round(x,2) for x in et]})
    print(f"{label:9s} K={K:<6}: base={base:5.2f} peak={et.max():5.2f} fold={fold:.3f} a*={astar:.2f} I_op={I_adhe(0,0.6,K):.2f}",flush=True)
json.dump({'alphas':list(map(float,alphas)),'rows':rows}, open(f'/home/claude/unified_data/adhe_{GSMM}.json','w'))
print("saved adhe_%s.json"%GSMM)
