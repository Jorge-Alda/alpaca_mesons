"""Plot finite conditional profile examples; no probability-band interpretation.
Run .venv/bin/python scripts/plot_odd_low_energy_profiles.py.
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def main():
    data=json.loads((ROOT/'data/vector_sector/odd_low_energy_profiles.json').read_text())
    assert data['complete_256_pair_run']
    controls=('primary','without_broad_rho','without_constrained_Kstar','narrow_nonconstrained_summary_control')
    titles=('Primary','Exclude broad rho','Exclude constrained K-star','Exclude both')
    fig,axes=plt.subplots(2,2,figsize=(9,6.5),sharex=True,sharey=True,layout='constrained')
    examples=[]
    for ax,control,title in zip(axes.flat,controls,titles):
        candidates=[r for r in data['results'] if r['control']==control and r['profiles'] and r['profile_reference_resolved']]
        ax.set_title(title);ax.axhline(1,color='gray',lw=.7);ax.axhline(3.84,color='gray',lw=.7,ls=':')
        t=np.linspace(-4,4,100);ax.plot(t,t*t,color='black',lw=1,ls='--',label='Local quadratic')
        if candidates:
            row=min(candidates,key=lambda r:r['best']['chi2']);examples.append(dict(control=control,scenario_index=row['scenario_index']))
            for i,profile in enumerate(row['profiles']):
                points=sorted(profile['points'],key=lambda p:p['standardized_target'])
                x=[p['standardized_target'] for p in points]
                y=[p['delta_chi2'] if p['valid'] else np.nan for p in points]
                ax.plot(x,y,'.-',lw=1,label=f'Physical contrast {i+1}')
            ax.text(.03,.96,f'Scenario {row["scenario_index"]}; achieved chi²={row["best"]["chi2"]:.3g}',transform=ax.transAxes,va='top',fontsize=8)
        else:ax.text(.5,.5,'No resolved profile example',transform=ax.transAxes,ha='center')
        ax.set_xlim(-4.5,4.5);ax.set_ylim(-.1,10);ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=8,loc='upper right')
    fig.supxlabel('Fixed physical-amplitude contrast / reference local standard deviation')
    fig.supylabel('Profile Δχ²')
    fig.suptitle('Finite conditional profiles; illustrative achieved minima',fontsize=11)
    for extension in ('png','pdf'):
        target=ROOT/f'data/vector_sector/odd_low_energy_profiles.{extension}'
        fig.savefig(target,dpi=170)
        print(target.relative_to(ROOT))
    (ROOT/'data/vector_sector/odd_low_energy_profile_plot_examples.json').write_text(json.dumps(examples,indent=2)+'\n')
    plt.close(fig)

if __name__=='__main__':main()
