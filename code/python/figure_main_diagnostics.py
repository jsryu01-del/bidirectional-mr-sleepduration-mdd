"""Main-analysis diagnostic figures (one per direction; panels a-d), as in Main_normal_liability_4_Figures:
  a: exposure-increasing-allele scatter with IVW, MR-Egger, weighted median and (if outliers) MR-PRESSO slopes
  b: funnel plot of single-SNP Wald ratios vs 1/SE
  c: single-SNP Wald ratios with 95% CI, followed by MR-Egger, weighted median and IVW estimates
  d: IVW leave-one-out estimates
Inputs : results/runs/published_seed/<analysis>/ and results/figure_data/ (written by code/R/01 and 02)
Outputs: figures/<analysis>_diagnostics.png (300 dpi) and .pdf
Usage  : python code/python/figure_main_diagnostics.py
Plotting code is taken unchanged from the published figure script; only input paths/column names were adapted.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[2]/'.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json

PKG=Path(__file__).resolve().parents[2]
ROOT=PKG/'results/runs/published_seed'
FIGDATA=PKG/'results/figure_data'
OUT=PKG/'figures'
OUT.mkdir(exist_ok=True)
FONT=12
LW=1.8
BOTTOM_PLOT_IN=4.2
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':FONT,'axes.labelsize':FONT,
 'xtick.labelsize':FONT,'ytick.labelsize':FONT,'legend.fontsize':FONT,'legend.title_fontsize':FONT,
 'axes.linewidth':1.2,'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
COLORS={'IVW':'#1B9E77','MR-Egger':'#D95F02','MR-PRESSO':'#B8860B','Weighted median':'#7570B3'}
FUNNEL_COLORS={'IVW':'#A6CEE3','MR-Egger':'#1F78B4','Weighted median':'#B2DF8A'}
FUNNEL_STYLE={'IVW':'-','MR-Egger':'--','Weighted median':'-.'}
STYLE={'IVW':'-','MR-Egger':(0,(6,4)),'MR-PRESSO':(0,(6,4)),'Weighted median':(0,(5,3,1,3))}
NAMES={'short':'short sleep duration','long':'long sleep duration','MDD':'major depressive disorder'}
def cap(t): return t[0].upper()+t[1:]  # sentence-initial only
DIRECTIONS=['short_to_MDD','long_to_MDD','MDD_to_short','MDD_to_long']
manifest=[]

def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=1.2,length=5,pad=3)
    ax.grid(False)

def forest(ax, data, labels, n_individual, red_start):
    y=np.arange(len(data))[::-1]
    # Bold lines for sparse panels; scale down only when rows are too dense to separate.
    row_pt=(BOTTOM_PLOT_IN*72)/len(data)
    lw_i=min(LW,max(.85,.5*row_pt)); lw_r=min(2.4,max(1.2,.6*row_pt))
    ms_i=min(4.5,max(2.5,.9*row_pt)); ms_r=min(6,max(3.5,1.1*row_pt))
    for i,(_,r) in enumerate(data.iterrows()):
        col='#FF0000' if i>=red_start else '#444444'
        ax.plot([r['lower'],r['upper']],[y[i],y[i]],color=col,lw=(lw_r if i>=red_start else lw_i),solid_capstyle='butt',zorder=2)
        ax.plot(r['b'],y[i],'o',ms=(ms_r if i>=red_start else ms_i),color=('#FF0000' if i>=red_start else 'black'),zorder=3)
    ax.axvline(0,color='black',lw=1.2,ls=(0,(4,3)),zorder=1)
    ax.axhline(y[red_start]+.5,color='#999999',lw=1.0)
    ax.set_yticks(y,labels,fontsize=FONT)
    ax.tick_params(axis='y',length=4,width=1.2)
    ax.set_ylim(-.8,len(data)-.2)
    ax.margins(x=.06)

for direction in DIRECTIONS:
    exposure,outcome=direction.split('_to_')
    dat=pd.read_csv(ROOT/direction/'input_IVs.tsv',sep='\t')
    effects=pd.read_csv(ROOT/direction/'causal_effects.tsv',sep='\t').rename(columns={'beta':'beta_logOR','ci_lower':'ci_lower_logOR','ci_upper':'ci_upper_logOR'}).set_index('method')
    intercept=pd.read_csv(ROOT/direction/'egger_intercept.tsv',sep='\t').iloc[0]['intercept']
    tests=pd.read_csv(ROOT/direction/'PRESSO_tests.tsv',sep='\t').iloc[0]
    ss=pd.read_csv(FIGDATA/f'{direction}_single_SNP.tsv',sep='\t')
    loo=pd.read_csv(FIGDATA/f'{direction}_leave_one_out.tsv',sep='\t')
    n=len(dat)
    # Equal text sizes and equal panel heights; rsID overlap is explicitly allowed.
    top_h=4.7
    bottom_h=top_h
    height=.35+top_h+1.05+bottom_h+1.0
    width=14
    bottom_y=1.0/height
    bottom_plot_h=(bottom_h-.5)/height
    top_y=(1.0+bottom_h+1.05)/height
    top_plot_h=(top_h-.5)/height
    fig=plt.figure(figsize=(width,height))
    axes=[fig.add_axes([.105,top_y,.29,top_plot_h]),fig.add_axes([.585,top_y,.255,top_plot_h]),
          fig.add_axes([.15,bottom_y,.32,bottom_plot_h]),fig.add_axes([.62,bottom_y,.31,bottom_plot_h])]
    for ax in axes: clean(ax)
    a,b,c,d=axes
    # Reorient SNPs to exposure-increasing alleles, matching the Egger fit.
    sign=np.where(dat['beta.exposure']<0,-1,1)
    x=dat['beta.exposure'].to_numpy()*sign
    y=dat['beta.outcome'].to_numpy()*sign
    a.scatter(x,y,s=22,c='#444444',linewidths=0,zorder=3)
    span=x.max()-x.min(); xx=np.array([max(0,x.min()-.06*span),x.max()+.06*span])
    line_methods=['IVW','MR-Egger','MR-PRESSO','Weighted median']
    if int(tests['n_outliers'])==0: line_methods.remove('MR-PRESSO')
    for method in line_methods:
        source='MR-PRESSO outlier-corrected' if method=='MR-PRESSO' else method
        assert effects.loc[source,'status']=='estimated'
        slope=effects.loc[source,'beta_logOR']
        offset=intercept if method=='MR-Egger' else 0
        a.plot(xx,offset+slope*xx,label=method,color=COLORS[method],ls=STYLE[method],lw=LW)
    a.set_xlim(xx)
    a.set_xlabel(f'SNP effect on {NAMES[exposure]}\n(log odds)')
    a.set_ylabel(f'SNP effect on {NAMES[outcome]}\n(log odds)')
    a.legend(title='MR method',frameon=False,loc='center left',bbox_to_anchor=(1.01,.5),borderaxespad=0,handlelength=2.5)
    b.scatter(ss['b'],1/ss['se'],s=22,c='black',linewidths=0)
    for method in ['IVW','MR-Egger','Weighted median']:
        # Distinct line types so IVW and weighted median stay distinguishable when their estimates overlap.
        b.axvline(effects.loc[method,'beta_logOR'],color=FUNNEL_COLORS[method],ls=FUNNEL_STYLE[method],lw=LW,label=method)
    b.set_xlabel(r'$\beta_{IV}$ (log OR per unit increase in exposure log-odds)')
    b.set_ylabel(r'$1/SE_{IV}$')
    b.legend(title='MR method',frameon=False,loc='center left',bbox_to_anchor=(1.01,.5),borderaxespad=0,handlelength=2.5)
    # Use saved method-specific intervals (Egger t, IVW and weighted median normal).
    indiv=ss.sort_values('b',ascending=False).copy()
    aggregate=[]
    for method in ['MR-Egger','Weighted median','IVW']:
        v=effects.loc[method]
        aggregate.append({'SNP':'All - '+method,'b':v.beta_logOR,'lower':v.ci_lower_logOR,'upper':v.ci_upper_logOR})
    cdata=pd.concat([indiv,pd.DataFrame(aggregate)],ignore_index=True)
    forest(c,cdata,cdata['SNP'],n,n)
    c.set_xlabel(f'MR effect (log OR per unit increase in exposure log-odds)\n{cap(NAMES[exposure])} → {NAMES[outcome]}')
    ldata=pd.concat([loo[loo.SNP!='All'].sort_values('b',ascending=False),loo[loo.SNP=='All']],ignore_index=True)
    forest(d,ldata,ldata['SNP'],n,n)
    d.set_xlabel(f'IVW leave-one-out effect (log OR per unit increase in exposure log-odds)\n{cap(NAMES[exposure])} → {NAMES[outcome]}')
    for label,ax in zip('abcd',axes):
        pos=ax.get_position()
        fig.text(.025 if label in 'ac' else .525,pos.y1+.12/height,label,fontsize=20,fontweight='bold',va='bottom')
    # Verify actual font sizes, selected model curves, and complete SNP lists.
    for ax in axes:
        text=[ax.xaxis.label,ax.yaxis.label,*ax.get_xticklabels(),*ax.get_yticklabels()]
        if ax.get_legend(): text += [ax.get_legend().get_title(),*ax.get_legend().get_texts()]
        assert all(t.get_fontsize()==FONT for t in text)
    assert ('MR-PRESSO' in line_methods)==(int(tests['n_outliers'])>0)
    assert set(cdata.SNP[:n])==set(dat.SNP) and set(ldata.SNP[:n])==set(dat.SNP)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for ax in axes:
        text=[ax.xaxis.label,ax.yaxis.label,*ax.get_xticklabels(),*ax.get_yticklabels()]
        if ax.get_legend(): text += [ax.get_legend().get_title(),*ax.get_legend().get_texts()]
        for t in text:
            if not t.get_text(): continue
            if t in ax.get_xticklabels() and not min(ax.get_xlim()) <= t.get_position()[0] <= max(ax.get_xlim()): continue
            if t in ax.get_yticklabels() and not min(ax.get_ylim()) <= t.get_position()[1] <= max(ax.get_ylim()): continue
            bb=t.get_window_extent(renderer)
            assert bb.x0>=0 and bb.x1<=fig.bbox.width and bb.y0>=0 and bb.y1<=fig.bbox.height, (direction,t.get_text(),bb)
    assert np.allclose([ax.get_position().height for ax in axes], top_plot_h)
    # Preserve all rsIDs at 9 pt, allowing overlap as requested.
    pd.DataFrame({'SNP':dat.SNP,'beta_exposure_oriented':x,'beta_outcome_oriented':y}).to_csv(FIGDATA/f'{direction}_scatter_data.tsv',sep='\t',index=False)
    cdata.to_csv(FIGDATA/f'{direction}_forest_data.tsv',sep='\t',index=False)
    for ext in ['png','pdf']:
        fig.savefig(OUT/f'{direction}_diagnostics.{ext}',dpi=300)
    # Separate zoomed top/bottom previews permit inspection without reducing the 174-row labels.
    manifest.append({'direction':direction,'nsnp':n,'width_inches':width,'height_inches':height,'dpi':300,
      'text_pt':FONT,'panel_tag_pt':20,'equal_panel_heights':True,'rsID_overlap_allowed':True,'PRESSO_line_in_a':'MR-PRESSO' in line_methods,'all_SNP_labels':True,'status':'PASS'})
    plt.close(fig)
pd.DataFrame(manifest).to_csv(FIGDATA/'plot_validation.tsv',sep='\t',index=False)
(FIGDATA/'style.json').write_text(json.dumps({'scatter_colors':COLORS,'funnel_colors':FUNNEL_COLORS,'text_size_pt':FONT,'panel_tag_size_pt':20},indent=2))
print(pd.DataFrame(manifest).to_string(index=False))
