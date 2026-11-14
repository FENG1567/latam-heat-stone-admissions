"""Publication plots from the verified revision freeze; no effect estimation."""
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text

HERE=Path(__file__).resolve().parent
OLD=HERE
sys.path.insert(0,str(HERE))
from publication_common import FrozenPublication,count,num,percent,source_csv
spec=importlib.util.spec_from_file_location('base_figure_style',HERE/'render_repaired_figures.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
BLUE,RED,GRAY=base.BLUE,base.RED,base.NEUTRAL

def role(pub,name,**where):return pub.rows(pub.metadata['revision']['tables_by_role'][name],**where)
def one_role(pub,name,**where):return pub.one(pub.metadata['revision']['tables_by_role'][name],**where)
def primary(pub,country,endpoint):
    return one_role(pub,'marginal_sync',country=country,endpoint=endpoint,bandwidth_days='14',
                    inference='calendar_synchronous_plus_DL_tau_component')
def stats(r):return f"k={count(r['k'])}, I²={num(r['I2'],1)}%, direction >0={percent(r['concordance'],1,fraction=True)}"
def ci(r):return f"{num(r['RR'])} ({num(r['lo'])}–{num(r['hi'])})"
def save(fig,path):
    for ext in ['pdf','svg','png']:
        fig.savefig(path.with_suffix('.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def fig2(pub,out,sources):
    rows=[primary(pub,c,e) for c in ['BR','MX'] for e in ['y_K80','y_urinary']]
    source_csv(rows,sources/'figure2_same_joint_sample_marginals.csv')
    labels=['Brazil | K80','Brazil | urinary N20–N23','Mexico | K80','Mexico | urinary N20–N23']
    fig,(ax,stat)=plt.subplots(1,2,figsize=(8.8,4.5),gridspec_kw={'width_ratios':[2.4,2.3]},constrained_layout=True)
    y=np.arange(4)[::-1]
    for i,(r,yy) in enumerate(zip(rows,y)):
        point,lo,hi=map(float,[r['RR'],r['lo'],r['hi']]);color=BLUE if r['country']=='BR' else RED
        ax.errorbar(point,yy,xerr=[[point-lo],[hi-point]],fmt='s' if r['endpoint']=='y_K80' else 'o',
                    color=color,elinewidth=1.3,capsize=3,markersize=5.5)
    ax.axvline(1,color=GRAY,ls='--',lw=1)
    low=min(float(r['lo']) for r in rows);high=max(float(r['hi']) for r in rows)
    ax.set_xlim(low-.015,high+.015);ax.set_ylim(-.7,3.8)
    ax.set_yticks(y,labels);ax.grid(axis='x',color=base.GRID);ax.set_xlabel('RR | local P95 versus shared median')
    ax.set_title('a  Same primary joint-fit samples',loc='left',fontweight='bold',fontsize=10)
    stat.axis('off');stat.set_ylim(ax.get_ylim());stat.set_xlim(0,1)
    columns=[0,.49,.60,.75]
    for x,t in zip(columns,['RR (95% CI)','k','I² (%)','Direction >0']):stat.text(x,3.55,t,fontsize=8,fontweight='bold')
    for r,yy in zip(rows,y):
        for x,t in zip(columns,[ci(r),count(r['k']),num(r['I2'],1),percent(r['concordance'],1,fraction=True)]):
            stat.text(x,yy,t,fontsize=8,va='center')
    stat.text(0,-.5,'DL point; synchronized HAC14 + tau uncertainty',fontsize=7.3,color=GRAY)
    fig.suptitle('Diagnosis-specific associations under the shared reference',fontweight='bold',fontsize=12)
    save(fig,out/'fig2_forest')
    return ('## Figure 2\n\n**Diagnosis-specific associations under the shared reference.** '
      'Four marginal RRs for K80 and broad urinary N20–N23 admission records use the same successful joint-fit units within each country. '
      'Every contrast compares the unit-specific P95 of E3 with its shared median. Bars show 95% confidence intervals from synchronized calendar-day Bartlett HAC14 plus the modeled DL between-unit component, conditional on observed weights and heterogeneity. '
      'These pooled marginal RRs use endpoint-specific pooling weights and must not be divided to reconstruct the directly estimated RRR. '
      'Colors identify Brazil (blue) and Mexico (red); squares identify K80 and circles urinary diagnoses. '
      'I² and the positive-unit proportion refer to unit HAC14 contrasts, not a new synchronized heterogeneity fit. '
      +' '.join(f"{label}: RR {ci(r)}, {stats(r)}." for label,r in zip(labels,rows)))

def fig4(pub,out,sources):
    rows=role(pub,'curves',country='BR',curve_scope='support_restricted')
    rows=[r for r in rows if float(r['temperature'])<=30]
    source_csv(rows,sources/'figure4_DL_REML_supported_curves.csv')
    specs=[('y_K80','K80',BLUE),('y_urinary','Urinary N20–N23',RED),('urinary_vs_K80','Direct urinary / K80',GRAY)]
    fig=plt.figure(figsize=(8.8,6.1));gs=fig.add_gridspec(2,3,height_ratios=[4.5,1.1],left=.07,right=.98,top=.86,bottom=.10,wspace=.32,hspace=.34)
    axes=[fig.add_subplot(gs[0,i]) for i in range(3)];track=fig.add_subplot(gs[1,:])
    counts=None;summaries=[]
    for ax,(endpoint,label,color),letter in zip(axes,specs,'abc'):
        for method,linestyle,alpha in [('DL','-',.17),('REML_profile_HK','--',.07)]:
            r=sorted([z for z in rows if z['endpoint']==endpoint and z['method']==method],key=lambda z:float(z['temperature']))
            if len(r)<2:raise ValueError(f'Missing curve {endpoint}/{method}')
            x=np.array([float(z['temperature']) for z in r]);point=np.array([float(z['RR']) for z in r])
            lo=np.array([float(z['lo']) for z in r]);hi=np.array([float(z['hi']) for z in r])
            ax.fill_between(x,lo,hi,color=color,alpha=alpha,lw=0)
            ax.plot(x,point,color=color,ls=linestyle,lw=1.6)
            this=([z['temperature'] for z in r],[z['k'] for z in r])
            if counts is None:counts=this
            elif this!=counts:raise ValueError('Curve methods/endpoints differ in support membership counts')
            if method=='DL':summaries.append(f"{label} at 30 °C: {stats(r[-1])}")
        ax.axhline(1,color=GRAY,ls=':',lw=.8);ax.set_xlim(14,30);ax.set_title(letter+'  '+label,loc='left',fontsize=9,fontweight='bold')
        ax.set_ylabel('RR' if endpoint!='urinary_vs_K80' else 'Direct RRR');ax.grid(color=base.GRID,lw=.6)
    track.step([float(x) for x in counts[0]],[int(x) for x in counts[1]],where='mid',color=GRAY,lw=1.4)
    track.set_xlim(14,30);track.set_ylim(bottom=0);track.grid(axis='y',color=base.GRID)
    track.set_xlabel('Three-day mean temperature (°C)');track.set_ylabel('Contributing units, k');track.set_title('d  Common local-support count',loc='left',fontsize=9,fontweight='bold')
    fig.legend(handles=[Line2D([],[],color=GRAY,label='DL',lw=1.6),Line2D([],[],color=GRAY,ls='--',label='Profiled REML–HK',lw=1.6)],
               loc='upper center',bbox_to_anchor=(.5,.94),ncol=2,frameon=False)
    fig.suptitle('Locally supported curves and pooling-method sensitivity',fontweight='bold',fontsize=12,y=.99)
    save(fig,out/'fig4_curves')
    support=pub.one('support_32c',country='BR')
    return ('## Figure 4\n\n**Locally supported curves and pooling-method sensitivity.** '
      'a,b, K80 and urinary N20–N23 RRs versus each region’s shared median. c, Direct joint urinary-versus-K80 RRR. '
      'Solid lines and darker bands show DL estimates and pointwise 95% intervals; dashed lines and lighter bands show profiled REML–Hartung–Knapp estimates and pointwise 95% intervals. '
      'd, The common contributor count for both methods and all endpoints. A region contributes only inside its E3 P2–P98 with at least 30 eligible days within ±1 °C. '
      'Changing geography, weights and heterogeneity estimates affect curve shape. Large sampling-variance ratios generated stability warnings; no variance floor or convenient-unit exclusion was applied. '
      'The lines do not trace one fixed national population, and a crossing or bend is not a biological threshold. Display ends at 30 °C. '
      +'; '.join(summaries)+'. '
      +f"At 32 °C, {count(support['regions_with_at_least_1_E3_day_within_1C'])}/{count(support['regions_assessed'])} regions and {count(support['total_region_days_within_1C_of_32C'])} region-days occur within ±1 °C; only {count(support['regions_marked_strictly_supported_at_32C'])} region meets the strict curve rule. This unsupported tail is not shown. "
      'Full temperature-specific k, I², directional proportions and weight diagnostics accompany the source data.')

def render(frozen, output_root):
    output_root=Path(output_root)
    out=output_root/'outputs';sources=out/'source_data';supp=out/'supplementary_figures'
    for p in [out,sources,supp]:p.mkdir(parents=True,exist_ok=True)
    pub=FrozenPublication(frozen);base._style();legends=[]
    save_original=base._save_figure
    def neutral_figure1_save(fig,path):
        def replace(text):
            return (text.replace('Mexico preferred','Mexico primary')
                        .replace('2019 restored','2019 included')
                        .replace('preferred 2015–2022 cohort','primary 2015–2022 cohort')
                        .replace('Admissions, not unique people','Administrative admission records'))
        for ax in fig.axes:
            if any('preferred' in t.get_text() for t in ax.get_yticklabels()):
                ax.set_yticks(ax.get_yticks(),[replace(t.get_text()) for t in ax.get_yticklabels()])
        for item in fig.findobj(match=Text):item.set_text(replace(item.get_text()))
        timeline=fig.axes[0]
        year_labels=timeline.get_xticklabels()
        year_labels[-2].set_ha('right');year_labels[-1].set_ha('left')
        for ax in fig.axes[1:]:
            panel=ax.texts[0];panel.set_position((-0.14,1.02))
            ax.xaxis.set_label_coords(0.44,-0.13)
        return save_original(fig,path)
    base._save_figure=neutral_figure1_save
    try:
        _,legend=base._figure1(pub,out,sources)
    finally:
        base._save_figure=save_original
    legend=legend.replace('Counts are admissions, not unique people.', 'Counts are administrative admission records, not adjudicated unique episodes or people.')
    legends.append(legend)
    # Preserve exploratory MRT estimates as a separately numbered supplement.
    _,mrtlegend=base._figure2(pub,supp,sources)
    for ext in ['pdf','svg','png']:
        (supp/('fig2_forest.'+ext)).replace(supp/('figS1_mrt_exploratory.'+ext))
    mrtlegend=mrtlegend.replace('## Figure 2','## Figure S1',1)+'\n\nThese are secondary exploratory selected-MRT results. Intervals condition on the selected P10–P60 minimum; unit coefficient simulation is not national selection-aware inference. Optimized-reference positive-unit proportions are not independent fixed-reference directional replication.'
    legends.append(fig2(pub,out,sources))
    original_one=base._one
    def replacement(f,key,artifact,**where):
        if artifact=='figure3' and key in ['core_primary_contrast','mx_updated_primary']:
            return primary(f,where['country'],'urinary_vs_K80')
        if artifact=='figure3' and where.get('scenario')=='temporal_replication_2023_fixed_2015_2022_reference':
            selection=f.metadata['revision']['display_selections'].get('mexico2023_followup')
            if selection:return f.one(selection['table_key'],**selection['where'])
        return original_one(f,key,artifact,**where)
    base._one=replacement
    try:_,legend=base._figure3(pub,out,sources)
    finally:base._one=original_one
    legend=legend.replace('including the 2023 interval, which spans 1','including the updated 2023 interval')
    legend=legend.replace('Points and bars show DerSimonian–Laird estimates and 95% confidence intervals;', 'Points show DerSimonian–Laird estimates. Bars show the stated 95% confidence intervals;')
    legend+='\n\nCountry-panel intervals use synchronized HAC14 plus the modeled DL component conditional on pooling weights and heterogeneity. Archived grid and period-refit intervals retain their stated unit-HAC14 DL inference. The 2023 refit includes available next-year discharge follow-up and remains association replication, not prediction validation.'
    legends.append(legend);legends.append(fig4(pub,out,sources))
    (out/'figure_legends.md').write_text('\n\n'.join(legends)+'\n',encoding='utf-8')
    (supp/'supplementary_figure_legends.md').write_text(mrtlegend+'\n',encoding='utf-8')
    pub.save_claims(out/'figure_source_rows_revision.json')
    manifest={'main_figures':4,'main_tables_contract':4,'maximum_curve_temperature_C':30,
              'frozen_path':pub.path.name,'figure_sources':'source_data','backend':'Python matplotlib',
              'notes':'PDF/SVG editable text; PNG300dpi. All effect uncertainty read from verified frozen records; no estimation in plotting.'}
    (out/'figure_revision_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest))
    return manifest

