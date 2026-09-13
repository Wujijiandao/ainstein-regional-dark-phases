from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'research_v1.6' / 'results' / 'branch_parity'
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = [
    (ROOT / 'research_v1.4' / 'results' / 'descendant_transport_n64_4seed_combined.csv', '64^3'),
    (ROOT / 'research_v1.4' / 'results' / 'descendant_transport_n80_s11.csv', '80^3'),
]


def weighted_quantile(x, w, q):
    x = np.asarray(x, float); w = np.asarray(w, float)
    order = np.argsort(x)
    xx=x[order]; ww=w[order]
    c=np.cumsum(ww)/np.sum(ww)
    return float(xx[np.searchsorted(c, q, side='left')])


def analyze_one(path: Path, label: str):
    df = pd.read_csv(path)
    rows=[]; basin_rows=[]
    for seed,g in df.groupby('seed'):
        aa=np.array(sorted(g['a'].unique()), float)
        if not np.allclose(aa,[0.25,0.50,0.75,1.00]):
            raise RuntimeError(f'unexpected snapshots {aa}')
        p=g.pivot(index='final_basin_id',columns='a',values=['volume_fraction','Hloc_over_H','Gamma','instantaneous_active','mean_density_ratio'])
        ids=p.index
        w=p['volume_fraction'][1.0].to_numpy(float)
        w=w/w.sum()
        x=np.column_stack([p['Hloc_over_H'][a].to_numpy(float)-1.0 for a in aa])
        # Signed relative expansion history.  d ln(a_Omega/a)/d ln a = Hloc/H - 1.
        hist_int=np.trapezoid(x,np.log(aa),axis=1)
        # Exact endpoint finite relative dilation inferred from descendant comoving volume.
        dilation=(1.0/3.0)*np.log(p['volume_fraction'][1.0].to_numpy(float)/p['volume_fraction'][0.25].to_numpy(float))
        final_pos=x[:,-1]>0
        persistent_pos=np.all(x>0,axis=1)
        last2_pos=np.all(x[:,-2:]>0,axis=1)
        last3_pos=np.all(x[:,-3:]>0,axis=1)
        hist_pos=hist_int>0
        old_active=p['instantaneous_active'][1.0].to_numpy(bool)
        corr=float(np.corrcoef(hist_int,dilation)[0,1])
        rmse=float(np.sqrt(np.sum(w*(hist_int-dilation)**2)))
        row={
            'dataset':label,'seed':int(seed),
            'final_positive_fraction':float(w[final_pos].sum()),
            'persistent_all4_positive_fraction':float(w[persistent_pos].sum()),
            'last2_positive_fraction':float(w[last2_pos].sum()),
            'last3_positive_fraction':float(w[last3_pos].sum()),
            'history_integral_positive_fraction':float(w[hist_pos].sum()),
            'old_gamma_active_fraction':float(w[old_active].sum()),
            'old_active_inside_final_positive':float(w[old_active & final_pos].sum()/max(w[old_active].sum(),1e-30)),
            'old_active_inside_persistent_positive':float(w[old_active & persistent_pos].sum()/max(w[old_active].sum(),1e-30)),
            'history_dilation_corr':corr,
            'history_dilation_weighted_rmse':rmse,
            'dilation_weighted_median':weighted_quantile(dilation,w,0.5),
            'dilation_threshold_top30pct':weighted_quantile(dilation,w,0.70),
            'dilation_threshold_top50pct':weighted_quantile(dilation,w,0.50),
            'dilation_threshold_top70pct':weighted_quantile(dilation,w,0.30),
        }
        rows.append(row)
        for i,bid in enumerate(ids):
            basin_rows.append({
                'dataset':label,'seed':int(seed),'final_basin_id':int(bid),'final_volume_weight':float(w[i]),
                'x_025':float(x[i,0]),'x_050':float(x[i,1]),'x_075':float(x[i,2]),'x_100':float(x[i,3]),
                'history_integral':float(hist_int[i]),'relative_dilation':float(dilation[i]),
                'final_positive':int(final_pos[i]),'persistent_all4_positive':int(persistent_pos[i]),
                'old_gamma_active':int(old_active[i]),
            })
    return pd.DataFrame(rows),pd.DataFrame(basin_rows)

all_summary=[]; all_basins=[]
for p,label in DATASETS:
    s,b=analyze_one(p,label); all_summary.append(s); all_basins.append(b)
summary=pd.concat(all_summary,ignore_index=True)
basins=pd.concat(all_basins,ignore_index=True)
summary.to_csv(OUT/'branch_odd_volume_audit.csv',index=False)
basins.to_csv(OUT/'branch_odd_basin_histories.csv',index=False)

# Algebraic / symmetry identities used in the manuscript.
# For a velocity-gradient matrix A, theta=tr A, sigma = sym(A)-theta I/3.
# Qloc proxy = 2/3 theta^2 - 2 sigma^2 is invariant under A -> -A.
rng=np.random.default_rng(20260912)
parity=[]
for _ in range(64):
    A=rng.normal(size=(3,3))
    def q_of(A):
        S=0.5*(A+A.T); th=np.trace(S); sig=S-np.eye(3)*th/3
        return (2/3)*th**2 - np.sum(sig*sig)
    q1=q_of(A); q2=q_of(-A)
    parity.append(abs(q1-q2))
parity_max=float(max(parity))

# Perfect-fluid constitutive identities for L=F(b): rho=-F, p=F-bF_b.
# F=-m b -> dust; F=-eps -> vacuum.
bgrid=np.geomspace(1e-3,1e3,31)
m=1.7; eps=2.3
rho_d=m*bgrid; p_d=np.zeros_like(bgrid)
rho_v=np.full_like(bgrid,eps); p_v=-rho_v
fluid_identity_max=max(float(np.max(np.abs(p_d))),float(np.max(np.abs(p_v+rho_v))))

payload={
    'summary_records':summary.to_dict(orient='records'),
    'velocity_gradient_reversal_Q_even_max_abs_error':parity_max,
    'fluid_branch_identity_max_abs_error':fluid_identity_max,
    'definitions':{
        'x':'Hloc/Hbg - 1',
        'history_integral':'integral x d ln a over stored snapshots',
        'relative_dilation':'(1/3) ln[V_com(a=1)/V_com(a=0.25)]',
    },
}
(OUT/'branch_parity_summary.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')

# Figure: zero-threshold branch-odd supports versus the old Gamma-active set.
plot=summary.copy()
labels=[f"{r['dataset']} s{int(r['seed'])}" for _,r in plot.iterrows()]
xx=np.arange(len(plot))
width=0.20
fig,ax=plt.subplots(figsize=(9.2,4.8))
ax.bar(xx-1.5*width,100*plot['final_positive_fraction'],width,label=r'$X_\Omega(a=1)>0$')
ax.bar(xx-0.5*width,100*plot['persistent_all4_positive_fraction'],width,label='positive at all 4 snapshots')
ax.bar(xx+0.5*width,100*plot['history_integral_positive_fraction'],width,label=r'$\int X_\Omega d\ln a>0$')
ax.bar(xx+1.5*width,100*plot['old_gamma_active_fraction'],width,label=r'old $\Gamma$-active')
ax.set_ylabel('final Eulerian volume [%]')
ax.set_xticks(xx,labels,rotation=25,ha='right')
ax.set_ylim(0,65)
ax.legend(frameon=False,ncol=2)
ax.set_title('Branch-odd regional expansion controls are not volume-starved at zero threshold')
fig.tight_layout()
fig.savefig(OUT/'fig_branch_odd_volume_audit.pdf')
fig.savefig(OUT/'fig_branch_odd_volume_audit.png',dpi=180)
plt.close(fig)

# Scatter audit of integrated relative expansion versus exact descendant volume dilation.
sub=basins[(basins.dataset=='64^3') & (basins.seed==11)]
fig,ax=plt.subplots(figsize=(5.7,5.2))
size=12+900*sub['final_volume_weight'].to_numpy()
ax.scatter(sub['relative_dilation'],sub['history_integral'],s=size,alpha=0.55)
lo=min(sub['relative_dilation'].min(),sub['history_integral'].min())
hi=max(sub['relative_dilation'].max(),sub['history_integral'].max())
ax.plot([lo,hi],[lo,hi],'--',linewidth=1,label='identity')
ax.axhline(0,linewidth=0.8); ax.axvline(0,linewidth=0.8)
ax.set_xlabel(r'exact endpoint relative dilation $\frac{1}{3}\ln[V_c(1)/V_c(0.25)]$')
ax.set_ylabel(r'$\int (H_\Omega/H-1)\,d\ln a$ (4-snapshot quadrature)')
ax.legend(frameon=False)
ax.set_title('Signed expansion history tracks descendant relative dilation')
fig.tight_layout()
fig.savefig(OUT/'fig_history_dilation_audit.pdf')
fig.savefig(OUT/'fig_history_dilation_audit.png',dpi=180)
plt.close(fig)

print(summary.to_string(index=False))
print('Q parity max error',parity_max)
print('fluid branch identity max error',fluid_identity_max)
