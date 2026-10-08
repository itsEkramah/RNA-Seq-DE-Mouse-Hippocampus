"""Data-driven, consistently labeled figures; no synthetic expression values."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list

from scripts.data import MARKERS, SAMPLES

WT = '#2874A6'
KO = '#C65B28'
GRAY = '#B8C1C9'
INK = '#243746'
COLORS = [WT]*3 + [KO]*3


def save(fig, folder, name, caption):
    fig.text(.02, .012, caption, fontsize=8, color=INK, va='bottom')
    fig.savefig(folder / f'{name}.png', dpi=300, facecolor='white', bbox_inches='tight')
    fig.savefig(folder / f'{name}.pdf', facecolor='white', bbox_inches='tight',
                metadata={'Creator':'GSE116773 reproducible analysis', 'CreationDate':None, 'ModDate':None})
    plt.close(fig)


def title(fig, heading, subtitle):
    fig.suptitle(heading, fontsize=16, fontweight='bold', x=.02, ha='left', y=.98, color=INK)
    fig.text(.02, .895, subtitle, fontsize=10, color=INK)


def make_figures(folder: Path, tables: Path, meta, tech, library, corr, norm, vst, res,
                 diagnostic, coords, variance, distances):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
                         'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.labelcolor':INK,'text.color':INK,'axes.edgecolor':'#71808C',
                         'xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42,
                         'axes.axisbelow':True})

    # All tested genes appear; FDR and effect size are explicitly separate criteria.
    fig, axes = plt.subplots(1,2,figsize=(12,5.8))
    fig.subplots_adjust(left=.08,right=.98,bottom=.17,top=.81,wspace=.28)
    title(fig,'Differential expression after Notch2 deletion',
          'Public GRCm38 counts | 3 biological samples per condition | CKO relative to wild type')
    groups = [('FDR >= 0.05',GRAY),('Not assigned an FDR','#E3E8EC'),
              ('Lower in CKO',WT),('Higher in CKO',KO)]
    for status,color in groups:
        block=res[res.status.eq(status)]
        axes[0].scatter(block.baseMean,block.log2FoldChange,s=7,c=color,alpha=.65,rasterized=True,
                        linewidths=0,label=f'{status} (n={len(block):,})')
        block=block[block.padj.notna()]
        axes[1].scatter(block.log2FoldChange,-np.log10(block.padj.clip(lower=np.finfo(float).tiny)),
                        s=9,c=color,alpha=.75,rasterized=True,linewidths=0)
    axes[0].set(xscale='log',xlabel='Mean normalized gene count',ylabel='Log2 fold change (CKO / WT)',title='A  Expression abundance and effect size')
    axes[1].set(xlabel='Log2 fold change (CKO / WT)',ylabel='-log10(BH-adjusted P-value)',title='B  Genome-wide statistical evidence')
    for y in [-1,1]: axes[0].axhline(y,color='#80909A',ls=':',lw=.8)
    axes[0].axhline(0,color=INK,lw=.7)
    axes[1].axhline(-np.log10(.05),color=INK,ls='--',lw=.8)
    for x in [-1,1]: axes[1].axvline(x,color='#80909A',ls=':',lw=.8)
    for sign in [-1,1]:
        top=res[res.log2FoldChange.mul(sign).gt(0) & res.padj.lt(.05)].nsmallest(2,'padj')
        for i,(gene,r) in enumerate(top.iterrows()):
            axes[1].annotate(gene,(r.log2FoldChange,-np.log10(max(r.padj,np.finfo(float).tiny))),
                             xytext=(sign*(12+i*9),12+i*18),textcoords='offset points',
                             fontsize=8,ha='left' if sign>0 else 'right',
                             arrowprops={'arrowstyle':'-','color':'#7C8993','lw':.6})
    axes[0].legend(fontsize=7,loc='best',frameon=False,markerscale=2)
    axes[1].margins(x=.18,y=.2)
    save(fig,folder,'differential_expression','Unshrunk model estimates. Color: FDR < 0.05; dotted lines: |log2 fold change| = 1.\nMissing adjusted P-values are omitted from panel B. Effect-size cutoffs are descriptive, not threshold tests.')

    fig,axes=plt.subplots(1,2,figsize=(11.5,5.6))
    fig.subplots_adjust(left=.09,right=.93,bottom=.18,top=.81,wspace=.42)
    title(fig,'Variation between biological samples','Blind variance-stabilizing transformation; no samples removed')
    pca_offsets={'WT1':(6,6),'WT2':(-30,-13),'WT3':(6,6),
                 'KO1':(6,16),'KO2':(6,6),'KO3':(6,-17)}
    for i,sample in enumerate(SAMPLES):
        axes[0].scatter(*coords.loc[sample],s=80,c=COLORS[i],marker='o' if i<3 else '^',edgecolor='white',lw=.7)
        axes[0].annotate(sample,coords.loc[sample],xytext=pca_offsets[sample],textcoords='offset points',fontsize=9)
    axes[0].set(xlabel=f'PC1 ({variance[0]:.1%} of variance)',ylabel=f'PC2 ({variance[1]:.1%} of variance)',title='A  PCA: 500 most variable genes')
    axes[0].margins(.22); axes[0].grid(alpha=.15)
    order=leaves_list(linkage(vst.T,method='average',metric='euclidean'))
    dist=distances.iloc[order,order]
    im=axes[1].imshow(dist,cmap='Blues',vmin=0)
    axes[1].set(xticks=range(6),yticks=range(6),xticklabels=dist.columns,yticklabels=dist.index,title='B  Sample distances: all retained genes')
    for i in range(6):
        for j in range(6):
            val=dist.iloc[i,j]
            axes[1].text(j,i,f'{val:.0f}',ha='center',va='center',fontsize=8,
                         color='white' if val>dist.to_numpy().max()*.6 else INK)
    fig.colorbar(im,ax=axes[1],fraction=.046,pad=.04,label='Euclidean distance')
    save(fig,folder,'sample_structure','WT = wild type; KO = Notch2 conditional knockout. Each point represents one biological RNA sample (pooled animals).\nPCA selection uses variance, not differential-expression significance. Distance order uses average-linkage clustering.')

    selected=res[res.padj.lt(.05)].sort_values(['padj','pvalue'],kind='stable').head(24).index
    if len(selected):
        centered=vst.loc[selected].sub(vst.loc[selected].mean(axis=1),axis=0)
        order=leaves_list(linkage(centered,method='average')) if len(selected)>1 else [0]
        centered=centered.iloc[order]
        centered.to_csv(tables/'heatmap_row_centered_vst.csv',float_format='%.12g',lineterminator='\n')
        fig,ax=plt.subplots(figsize=(8,9))
        fig.subplots_adjust(left=.26,right=.82,bottom=.13,top=.84)
        title(fig,'Expression patterns of the top DE genes',f'Top {len(selected)} genes ranked by adjusted P-value (FDR < 0.05)')
        vmax=float(np.abs(centered.to_numpy()).max())
        im=ax.imshow(centered,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax))
        ax.set(yticks=range(len(centered)),yticklabels=centered.index,xticks=range(6),xticklabels=SAMPLES)
        for label,color in zip(ax.get_xticklabels(),COLORS): label.set_color(color)
        ax.axvline(2.5,color='white',lw=2)
        fig.colorbar(im,ax=ax,fraction=.045,pad=.05,label='VST expression minus gene mean')
        save(fig,folder,'top_gene_heatmap','Rows are centered VST values, not z-scores. Columns retain biological sample order.\nRows use average-linkage clustering. Selection by DE makes this a descriptive display, not independent validation.')

    fig,axes=plt.subplots(2,4,figsize=(12,8))
    fig.subplots_adjust(left=.07,right=.98,bottom=.14,top=.81,wspace=.42,hspace=.5)
    title(fig,'Candidate genes in the Notch2 study','A literature-motivated panel shown in full, regardless of statistical significance')
    marker_long=[]
    for gene,ax in zip(MARKERS,axes.flat):
        if gene not in norm.index:
            ax.text(.5,.5,'Below expression filter',ha='center',transform=ax.transAxes); ax.set_title(gene); continue
        for group,ids,color in [(0,SAMPLES[:3],WT),(1,SAMPLES[3:],KO)]:
            values=np.log2(norm.loc[gene,ids].to_numpy()+1)
            ax.scatter(group+np.array([-.13,0,.13]),values,c=color,s=42,edgecolors='white',lw=.6,zorder=3)
            ax.plot([group-.22,group+.22],[np.median(values)]*2,color=color,lw=1.5)
            for sample,val,logval in zip(ids,norm.loc[gene,ids],values):
                marker_long.append({'gene':gene,'sample':sample,'condition':'WildType' if group==0 else 'Notch2CKO',
                                    'normalized_count':val,'log2_normalized_count_plus_1':logval})
        q=res.loc[gene,'padj']
        ax.set_title(f'{gene}\nFDR = {q:.3g}' if pd.notna(q) else f'{gene}\nFDR unavailable',fontsize=11)
        ax.set(xticks=[0,1],xticklabels=['WT','CKO'],xlim=(-.5,1.5),ylabel='log2(normalized count + 1)')
        ax.grid(axis='y',alpha=.15)
    pd.DataFrame(marker_long).to_csv(tables/'candidate_gene_counts.csv',index=False,float_format='%.12g',lineterminator='\n')
    save(fig,folder,'candidate_gene_counts','Three biological samples per group; horizontal bars mark medians, and every point is shown. Each panel has its own y-axis range.\nPanel choice follows the source study and its analysis code. Expression associations alone do not establish direct regulation.')

    fig,ax=plt.subplots(figsize=(9,6.8))
    fig.subplots_adjust(left=.16,right=.79,bottom=.2,top=.80)
    title(fig,'Effect sizes and uncertainty for candidate genes','Notch2 CKO relative to wild type; same eight-gene panel')
    for i,gene in enumerate(MARKERS):
        if gene not in res.index: continue
        row=res.loc[gene]; color=KO if row.log2FoldChange>0 else WT
        ax.errorbar(row.log2FoldChange,i,xerr=1.959963984540054*row.lfcSE,fmt='o',color=color,capsize=3,ms=6)
        ax.text(1.04,i, f'{row.padj:.3g}' if pd.notna(row.padj) else 'NA',transform=ax.get_yaxis_transform(),fontsize=9,va='center')
    ax.text(1.04,1.03,'Genome-wide FDR',transform=ax.transAxes,fontsize=9)
    ax.axvline(0,color=INK,lw=.8)
    ax.set(yticks=range(8),yticklabels=MARKERS,xlabel='Log2 fold change (CKO / WT)',ylim=(7.6,-.6))
    ax.grid(axis='x',alpha=.15)
    save(fig,folder,'candidate_gene_effects','Points: unshrunk maximum-likelihood effects. Bars: approximate 95% Wald confidence intervals, not multiplicity-adjusted.\nFDR values use all eligible genes, not only this panel. With n = 3 per group, uncertainty and sample variability remain important.')

    fig,axes=plt.subplots(1,2,figsize=(11.5,5.8))
    fig.subplots_adjust(left=.08,right=.97,bottom=.18,top=.81,wspace=.3)
    title(fig,'Model diagnostics','Negative-binomial differential-expression model fitted to biological samples')
    axes[0].scatter(res.baseMean,diagnostic.loc[res.index,'genewise_dispersions'],c=GRAY,s=4,alpha=.5,rasterized=True,label='Gene-wise')
    axes[0].scatter(res.baseMean,diagnostic.loc[res.index,'dispersions'],c=WT,s=4,alpha=.45,rasterized=True,label='Final')
    sort=res.baseMean.sort_values().index
    axes[0].plot(res.loc[sort,'baseMean'],diagnostic.loc[sort,'fitted_dispersions'],c=KO,lw=1.5,label='Fitted trend')
    axes[0].set(xscale='log',yscale='log',xlabel='Mean normalized gene count',ylabel='Dispersion',title='A  Dispersion estimates')
    axes[0].legend(frameon=False,fontsize=8,markerscale=3)
    axes[1].hist(res.pvalue.dropna(),bins=np.linspace(0,1,41),color=WT,edgecolor='white',lw=.4)
    axes[1].set(xlabel='Raw Wald P-value',ylabel='Number of genes',title='B  Distribution of tested P-values',xlim=(0,1))
    save(fig,folder,'model_diagnostics','Missing P-values remain missing and are omitted from the histogram. The first bin contains P-values from 0 to 0.025.\nThese diagnostics assess model behavior; they do not verify read alignment or replace sample-level quality control.')

    ordered=[f'{s}_rep{lib}_{lane}' for s in SAMPLES for lib in [1,2] for lane in [1,2]]
    fig,(ax0,ax1)=plt.subplots(2,1,figsize=(11,12),gridspec_kw={'height_ratios':[1,3]})
    fig.subplots_adjust(left=.15,right=.9,bottom=.14,top=.88,hspace=.34)
    title(fig,'Technical files and their biological origin','Two library preparations x two sequencing lanes per biological RNA sample')
    x=np.arange(24); colors=[WT]*12+[KO]*12
    ax0.bar(x,library.loc[ordered,'assigned_gene_counts']/1e6,color=colors,width=.8)
    ax0.set(xticks=np.arange(1.5,24,4),xticklabels=SAMPLES,ylabel='Assigned gene counts (millions)',title='A  Sequencing contribution from each technical file')
    im=ax1.imshow(corr.loc[ordered,ordered],cmap='viridis',vmin=0,vmax=1)
    ax1.set(xticks=x,yticks=x,xticklabels=ordered,yticklabels=ordered,title='B  Technical-file expression correlations')
    ax1.tick_params(axis='x',labelrotation=90,labelsize=7); ax1.tick_params(axis='y',labelsize=7)
    for boundary in np.arange(3.5,23,4):
        ax1.axhline(boundary,c='white',lw=1); ax1.axvline(boundary,c='white',lw=1)
    fig.colorbar(im,ax=ax1,fraction=.035,pad=.03,label='Pearson r on log2(CPM + 1)')
    save(fig,folder,'technical_qc','CPM is used only for this descriptive technical comparison. Differential expression uses raw integer counts summed within each biological sample.\nBlue: wild type; orange: Notch2 CKO. Correlation scale spans 0 to 1; technical files are not independent biological replicates.')
