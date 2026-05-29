"""Unified dFBA on BOTH uploaded models. Optimized: configure each chassis-GSMM once, reuse."""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, cobra, time

class P:
    H2O2_0=5.0; NO2_0=2.5; NO3_0=1.0; ONOO_0=0.5
    k_H2O2=0.277; k_NO2=0.150; k_NO3=0.020; k_ONOO=1.500
    w_H2O2=1.0; w_NO2=0.4; w_NO3=0.05; w_ONOO=5.0
    kd_max=4.0; S_half_d=14.0; n_d=4; beta_ATPM=0.10
    r_base=0.5; r_max=10.0; Ka_r=1.8; n_a=2; Ki_r=13.0; n_i=4
    PFL_inh_K=0.8; PFL_inh_n=3; gamma_NADH=0.02
    ADH_eng=-40.0; dt=0.25; T_batch=10.0
    glc_0=55.0; X_0=0.05; glc_uptake_max=10.0

def stress_S(t,a):
    return (P.w_H2O2*a*P.H2O2_0*np.exp(-P.k_H2O2*t) + P.w_NO2*a*P.NO2_0*np.exp(-P.k_NO2*t)
            + P.w_NO3*a*P.NO3_0*np.exp(-P.k_NO3*t) + P.w_ONOO*a*P.ONOO_0*np.exp(-P.k_ONOO*t))
def r_func(S):
    return P.r_base + (P.r_max-P.r_base)*((S**P.n_a)/(P.Ka_r**P.n_a+S**P.n_a))*(1.0/(1.0+(S/P.Ki_r)**P.n_i))
def pfl_inhibition(S): return 1.0/(1.0+(S/P.PFL_inh_K)**P.PFL_inh_n)
def kd_func(S): return P.kd_max*(S**P.n_d)/(P.S_half_d**P.n_d+S**P.n_d)
def atpm_burden(S): return P.beta_ATPM*S + P.gamma_NADH*S

def load_model(path):
    m = cobra.io.read_sbml_model(path)
    m.reactions.EX_o2_e.lower_bound = 0.0
    m.reactions.EX_glc__D_e.lower_bound = -10.0
    return m

def biomass_id(m):
    for r in m.reactions:
        if 'BIOMASS' in r.id.upper() and 'core' in r.id.lower(): return r.id
    for r in m.reactions:
        if 'BIOMASS' in r.id.upper(): return r.id

def make_configured(model, chassis):
    """Configure a chassis ONCE; returns a model with the coupling constraint added."""
    m = model.copy()
    m.objective = biomass_id(m)
    m.reactions.EX_o2_e.lower_bound = 0.0
    m.reactions.EX_glc__D_e.lower_bound = -P.glc_uptake_max
    if chassis in ('KO11','LY160'):
        for rid in ['FRD2','FRD3']:
            if rid in m.reactions: m.reactions.get_by_id(rid).bounds=(0,0)
        if 'LDH_D' in m.reactions: m.reactions.LDH_D.bounds=(0,0)
        m.reactions.ALCD2x.lower_bound = P.ADH_eng
    if chassis == 'LY160':
        if 'ACKr' in m.reactions: m.reactions.ACKr.bounds=(0,0)
    cc = m.problem.Constraint(m.reactions.EX_etoh_e.flux_expression
                              - P.r_base*m.reactions.EX_for_e.flux_expression,
                              lb=0, ub=None, name='redox_coupling')
    m.add_cons_vars([cc])
    return m

def set_coupling(m, rval):
    c = m.constraints['redox_coupling']
    c.set_linear_coefficients({
        m.reactions.EX_etoh_e.forward_variable: 1.0, m.reactions.EX_etoh_e.reverse_variable:-1.0,
        m.reactions.EX_for_e.forward_variable:-rval, m.reactions.EX_for_e.reverse_variable: rval})

def run_dfba_cfg(m, alpha, scenario='full'):
    """Run dFBA on a PRE-CONFIGURED model m (reused). Resets bounds each call."""
    dt=P.dt; n=int(P.T_batch/dt)
    pfl0 = m.reactions.PFL.upper_bound
    # ATPM baseline: capture the model's native NGAM once
    atpm0 = getattr(m, '_atpm0', None)
    if atpm0 is None:
        atpm0 = m.reactions.ATPM.lower_bound; m._atpm0 = atpm0
    alcd_lb0 = m.reactions.ALCD2x.lower_bound
    if scenario=='no_alcd': m.reactions.ALCD2x.bounds=(0,0)
    X=P.X_0; glc=P.glc_0; etoh=for_=ac=0.0; pk=0.0
    for i in range(n):
        t=i*dt; S=stress_S(t,alpha); pk=max(pk,S)
        if scenario=='baseline': rv=P.r_base; pf=1.0
        elif scenario=='coupling': rv=r_func(S); pf=1.0
        elif scenario=='pfl': rv=P.r_base; pf=pfl_inhibition(S)
        else: rv=r_func(S); pf=pfl_inhibition(S)
        set_coupling(m, rv)
        m.reactions.PFL.upper_bound = pfl0*pf
        m.reactions.ATPM.lower_bound = atpm0 + atpm_burden(S)
        m.reactions.EX_glc__D_e.lower_bound = 0.0 if glc<0.1 else -P.glc_uptake_max
        try:
            sol=m.optimize()
            if sol.status!='optimal': X=X*np.exp(-kd_func(S)*dt); continue
            X=X*np.exp((sol.objective_value-kd_func(S))*dt)
            glc=max(0,glc+sol.fluxes['EX_glc__D_e']*X*dt)
            etoh+=sol.fluxes['EX_etoh_e']*X*dt; for_+=sol.fluxes['EX_for_e']*X*dt; ac+=sol.fluxes['EX_ac_e']*X*dt
        except Exception: pass
    # restore
    m.reactions.PFL.upper_bound=pfl0; m.reactions.ATPM.lower_bound=atpm0
    m.reactions.ALCD2x.lower_bound=alcd_lb0
    gu=P.glc_0-glc
    return {'EtOH':etoh,'For':for_,'Ac':ac,'X':X,'peak_S':pk,'yield':etoh/gu if gu>0 else 0,'glc':glc}
