#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import gc 
import os


sc.settings.verbosity = 3             
sc.settings.set_figure_params(dpi=100,figsize=(4,4))


# # BulkRNA on Rat

# In[2]:


p4vsp60 = pd.read_csv('../../../../../bulkRNA/DE/DE_markers/P4vsP60.txt',header=0,delimiter='\t')
p10vsp60 = pd.read_csv('../../../../../bulkRNA/DE/DE_markers/P10vsP60.txt',header=0,delimiter='\t')
p4vsp10 = pd.read_csv('../../../../../bulkRNA/DE/DE_markers/P10vsP60.txt',header=0,delimiter='\t')


# ## Volcano plots with original rat genes

# In[3]:


from bioinfokit import analys, visuz


# ### P4 vs P60

# In[4]:


p4vsp60 = p4vsp60[~p4vsp60.pvalue.isnull()]
p4vsp60['GeneNames'] = p4vsp60.index.tolist()
p4vsp60_sig  = p4vsp60[~p4vsp60.padj.isnull()]
p4vsp60_sig = p4vsp60_sig[(p4vsp60_sig.padj<0.05)&(abs(p4vsp60_sig.log2FoldChange)>2)]
p4vsp60_sig = p4vsp60_sig.sort_values("padj" ,ascending=True)
visuz.GeneExpression.volcano(df=p4vsp60, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='p4vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.4,1),
                             r=300,show=True,genenames =tuple(p4vsp60_sig['GeneNames'][:50].tolist()))


# ### P4 vs P10

# In[5]:


p4vsp10 = p4vsp10[~p4vsp10.pvalue.isnull()]
p4vsp10['GeneNames'] = p4vsp10.index.tolist()
p4vsp10_sig  = p4vsp10[~p4vsp10.padj.isnull()]
p4vsp10_sig = p4vsp10_sig[(p4vsp10_sig.padj<0.05)&(abs(p4vsp10_sig.log2FoldChange)>2)]
p4vsp10_sig = p4vsp10_sig.sort_values("padj" ,ascending=True)
visuz.GeneExpression.volcano(df=p4vsp10, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='p4vsp10',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.4,1),
                             r=300,show=True,genenames =tuple(p4vsp10_sig['GeneNames'][:50].tolist()))


# ### P10 vs P60

# In[6]:


p10vsp60 = p10vsp60[~p10vsp60.pvalue.isnull()]
p10vsp60['GeneNames'] = p10vsp60.index.tolist()
p10vsp60_sig  = p10vsp60[~p10vsp60.padj.isnull()]
p10vsp60_sig = p10vsp60_sig[(p10vsp60_sig.padj<0.05)&(abs(p10vsp60_sig.log2FoldChange)>2)]
p10vsp60_sig = p10vsp60_sig.sort_values("padj" ,ascending=True)
visuz.GeneExpression.volcano(df=p10vsp60, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='p10vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.4,1),
                             r=300,show=True,genenames =tuple(p10vsp60_sig['GeneNames'][:50].tolist()))


# # Convert Rat genes to mouse genes

# In[7]:


ConvertedGenelist = pd.read_table('RatGeneConvertedMouseGenes_list.txt',index_col=0)
ConvertedGenelist_edit = ConvertedGenelist.drop_duplicates('RGD.symbol', keep='last')
ConvertedGenelist_edit['GeneNames'] = ConvertedGenelist_edit['RGD.symbol']


# In[8]:


#P4 vs P60
p4vsp60_combine =  pd.merge(p4vsp60, ConvertedGenelist_edit, how='inner')
p4vsp60_combine = p4vsp60_combine[~p4vsp60_combine['MGI.symbol'].isnull()]
#P10 vs P60
p10vsp60_combine =  pd.merge(p10vsp60, ConvertedGenelist_edit, how='inner')
p10vsp60_combine = p10vsp60_combine[~p10vsp60_combine['MGI.symbol'].isnull()]
#P4 vs P10
p4vsp10_combine =  pd.merge(p4vsp10, ConvertedGenelist_edit, how='inner')
p4vsp10_combine = p4vsp10_combine[~p4vsp10_combine['MGI.symbol'].isnull()]


# # Mesenchymal

# In[9]:


mesen_DE = sc.read_h5ad('../results/Mesenchymal_DE.h5ad')


# In[10]:


sc.tl.rank_genes_groups(mesen_DE, groupby='Age', groups=['P3','P60'],reference='P60',method = 'wilcoxon', key_added = "wilcoxon_P3vsP60",use_raw=False)
sc.tl.rank_genes_groups(mesen_DE, groupby='Age', groups=['P3','P60'],reference='P3',method = 'wilcoxon', key_added = "wilcoxon_P60vsP3",use_raw=False)

sc.tl.rank_genes_groups(mesen_DE, groupby='Age', groups=['P11','P60'],reference='P60',method = 'wilcoxon', key_added = "wilcoxon_P11vsP60",use_raw=False)
sc.tl.rank_genes_groups(mesen_DE, groupby='Age', groups=['P11','P60'],reference='P11',method = 'wilcoxon', key_added = "wilcoxon_P60vsP11",use_raw=False)

sc.tl.rank_genes_groups(mesen_DE, groupby='Age', groups=['P3','P11'],reference='P11',method = 'wilcoxon', key_added = "wilcoxon_P3vsP11",use_raw=False)
sc.tl.rank_genes_groups(mesen_DE, groupby='Age', groups=['P3','P11'],reference='P3',method = 'wilcoxon', key_added = "wilcoxon_P11vsP3",use_raw=False)

sc.tl.rank_genes_groups(mesen_DE, groupby='Age', method = 'wilcoxon', key_added = "wilcoxon_Age",use_raw=False)


# ## P3 vs P60

# In[11]:


p3vsp60 = sc.get.rank_genes_groups_df(mesen_DE,group='P3',key='wilcoxon_P3vsP60')
p3vsp60 =p3vsp60.set_index(keys='names')
p3vsp60_sig = p3vsp60[(p3vsp60['pvals_adj']<0.05)& (abs(p3vsp60['logfoldchanges'])>2)]


# ### DEs identified in both Bulk and scRNA
# P4 vs P60

# In[12]:


p4vsp60_combine_sig = p4vsp60_combine[(abs(p4vsp60_combine['log2FoldChange'])>2)&(p4vsp60_combine['padj']<0.05)]

P3P60_overlap = []
for gene in p4vsp60_combine_sig['MGI.symbol'].tolist():
    if gene in p3vsp60_sig.index.tolist():
        P3P60_overlap.append(gene)

print('The number of overlapped genes: ',len(P3P60_overlap))


# #### Visualize in scRNA 

# In[13]:


for i in np.arange(6,len(P3P60_overlap)+6,6):
    sc.pl.violin(mesen_DE,keys=P3P60_overlap[(i-6):i],ncols=6,wspace=0.1,groupby='Age',rotation=90,multi_panel=True,jitter=0)


# #### Visualize in BulkRNA

# In[14]:


visuz.GeneExpression.volcano(df=p4vsp60_combine, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='4vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.42,1),
                             r=300,show=True,genenames =tuple(P3P60_overlap))


# ## P11 vs P60

# In[15]:


p11vsp60 = sc.get.rank_genes_groups_df(mesen_DE,group='P11',key='wilcoxon_P11vsP60')
p11vsp60 =p11vsp60.set_index(keys='names')
p11vsp60_sig = p11vsp60[(p11vsp60['pvals_adj']<0.05)& (abs(p11vsp60['logfoldchanges'])>2)]


# ### Identify DEs identified in both Bulk and scRNA
# P10 vs P60

# In[16]:


p10vsp60_combine_sig = p10vsp60_combine[(abs(p10vsp60_combine['log2FoldChange'])>2)&(p10vsp60_combine['padj']<0.05)]

P11P60_overlap = []
for gene in p10vsp60_combine_sig['MGI.symbol'].tolist():
    if gene in p11vsp60_sig.index.tolist():
        P11P60_overlap.append(gene)

print('The number of overlapped genes: ',len(P11P60_overlap))


# #### Visualize in scRNA 

# In[17]:


for i in np.arange(6,len(P11P60_overlap)+6,6):
    sc.pl.violin(mesen_DE,keys=P11P60_overlap[(i-6):i],ncols=6,wspace=0.1,groupby='Age',rotation=90,multi_panel=True,jitter=0)


# #### Visualize in BulkRNA

# In[18]:


visuz.GeneExpression.volcano(df=p10vsp60_combine, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='p10vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.42,1),
                             r=300,show=True,genenames =tuple(P11P60_overlap))


# ## P3 vs P11

# In[19]:


p3vsp11 = sc.get.rank_genes_groups_df(mesen_DE,group='P3',key='wilcoxon_P3vsP11')
p3vsp11 = p3vsp11.set_index(keys='names')
p3vsp11_sig = p3vsp11[(p3vsp11['pvals_adj']<0.05)& (abs(p3vsp11['logfoldchanges'])>2)]


# ### Identify DEs identified in both Bulk and scRNA
# P4 vs P10

# In[20]:


p4vsp10_combine_sig = p4vsp10_combine[(abs(p4vsp10_combine['log2FoldChange'])>2)&(p4vsp10_combine['padj']<0.05)]

P3P11_overlap = []
for gene in p4vsp10_combine_sig['MGI.symbol'].tolist():
    if gene in p3vsp11_sig.index.tolist():
        P3P11_overlap.append(gene)

print('The number of overlapped genes: ',len(P3P11_overlap))


# No overlaped genes in P3 vs P11

# ## Overlap summary
# Order genes by up-regulation and down-regulation in P3 or P11

# In[21]:


#the total overlapped genes
overlap_mesen = P3P60_overlap+P11P60_overlap
overlap_mesen = list(set(overlap_mesen))
print('the total number of DEs identified in both data set:',len(overlap_mesen))

overlap_mesen_df = p3vsp60.loc[overlap_mesen]
overlap_mesen_df=overlap_mesen_df.sort_values('scores',ascending=False)

sc.pl.heatmap(mesen_DE,var_names=overlap_mesen_df.index.tolist(),groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(12,6))



# ## Unique to scRNA

# In[22]:


sc.set_figure_params(scanpy=True, fontsize=10)


# In[23]:


p3vsp60_unique = []
for gene in p3vsp60_sig.index.tolist():
    if gene not in overlap_mesen:
        p3vsp60_unique.append(gene)
print('The number of DEs only found in scRNA between P3 vs P60: ',len(p3vsp60_unique))

p11vsp60_unique = []
for gene in p11vsp60_sig.index.tolist():
    if gene not in overlap_mesen:
        p11vsp60_unique.append(gene)
print('The number of DEs only found in scRNA between P11 vs P60: ',len(p11vsp60_unique))

p3vsp11_unique = []
for gene in p3vsp11_sig.index.tolist():
    if gene not in overlap_mesen:
        p3vsp11_unique.append(gene)
print('The number of DEs only found in scRNA between P3 vs P11: ',len(p3vsp11_unique))

#the total unique DEs in scRNA
scrna_DEs_mesen = list(set(p3vsp60_unique+p11vsp60_unique+p3vsp11_unique))
print('The total number of DEs found in scRNA: ',len(scrna_DEs_mesen))
scrna_DEs_mesen_df = p3vsp60.loc[scrna_DEs_mesen]
scrna_DEs_mesen_df=scrna_DEs_mesen_df.sort_values('scores',ascending=False)

sc.pl.heatmap(mesen_DE,var_names=scrna_DEs_mesen_df.index.tolist(),groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(24,6))


# # Endothelial

# In[24]:


endoth_DE = sc.read_h5ad('../results/Endothelial_DE.h5ad')


# In[25]:


sc.tl.rank_genes_groups(endoth_DE, groupby='Age', groups=['P3','P60'],reference='P60',method = 'wilcoxon', key_added = "wilcoxon_P3vsP60",use_raw=False)
sc.tl.rank_genes_groups(endoth_DE, groupby='Age', groups=['P3','P60'],reference='P3',method = 'wilcoxon', key_added = "wilcoxon_P60vsP3",use_raw=False)

sc.tl.rank_genes_groups(endoth_DE, groupby='Age', groups=['P11','P60'],reference='P60',method = 'wilcoxon', key_added = "wilcoxon_P11vsP60",use_raw=False)
sc.tl.rank_genes_groups(endoth_DE, groupby='Age', groups=['P11','P60'],reference='P11',method = 'wilcoxon', key_added = "wilcoxon_P60vsP11",use_raw=False)

sc.tl.rank_genes_groups(endoth_DE, groupby='Age', groups=['P3','P11'],reference='P11',method = 'wilcoxon', key_added = "wilcoxon_P3vsP11",use_raw=False)
sc.tl.rank_genes_groups(endoth_DE, groupby='Age', groups=['P3','P11'],reference='P3',method = 'wilcoxon', key_added = "wilcoxon_P11vsP3",use_raw=False)

sc.tl.rank_genes_groups(endoth_DE, groupby='Age', method = 'wilcoxon', key_added = "wilcoxon_Age",use_raw=False)


# ## P3 vs P60

# In[26]:


p3vsp60 = sc.get.rank_genes_groups_df(endoth_DE,group='P3',key='wilcoxon_P3vsP60')
p3vsp60 =p3vsp60.set_index(keys='names')
p3vsp60_sig = p3vsp60[(p3vsp60['pvals_adj']<0.05)& (abs(p3vsp60['logfoldchanges'])>2)]


# ### DEs identified in both Bulk and scRNA
# P4 vs P60

# In[27]:


p4vsp60_combine_sig = p4vsp60_combine[(abs(p4vsp60_combine['log2FoldChange'])>2)&(p4vsp60_combine['padj']<0.05)]

P3P60_overlap = []
for gene in p4vsp60_combine_sig['MGI.symbol'].tolist():
    if gene in p3vsp60_sig.index.tolist():
        P3P60_overlap.append(gene)

print('The number of overlapped genes: ',len(P3P60_overlap))


# #### Visualize in scRNA 

# In[28]:


for i in np.arange(6,len(P3P60_overlap)+6,6):
    sc.pl.violin(endoth_DE,keys=P3P60_overlap[(i-6):i],ncols=6,wspace=0.1,groupby='Age',rotation=90,multi_panel=True,jitter=0)


# #### Visualize in BulkRNA

# In[29]:


visuz.GeneExpression.volcano(df=p4vsp60_combine, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='4vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.42,1),
                             r=300,show=True,genenames =tuple(P3P60_overlap))


# ## P11 vs P60

# In[30]:


p11vsp60 = sc.get.rank_genes_groups_df(endoth_DE,group='P11',key='wilcoxon_P11vsP60')
p11vsp60 =p11vsp60.set_index(keys='names')
p11vsp60_sig = p11vsp60[(p11vsp60['pvals_adj']<0.05)& (abs(p11vsp60['logfoldchanges'])>2)]


# ### Identify DEs identified in both Bulk and scRNA
# P10 vs P60

# In[31]:


p10vsp60_combine_sig = p10vsp60_combine[(abs(p10vsp60_combine['log2FoldChange'])>2)&(p10vsp60_combine['padj']<0.05)]

P11P60_overlap = []
for gene in p10vsp60_combine_sig['MGI.symbol'].tolist():
    if gene in p11vsp60_sig.index.tolist():
        P11P60_overlap.append(gene)

print('The number of overlapped genes: ',len(P11P60_overlap))


# #### Visualize in scRNA 

# In[32]:


for i in np.arange(6,len(P11P60_overlap)+6,6):
    sc.pl.violin(endoth_DE,keys=P11P60_overlap[(i-6):i],ncols=6,wspace=0.1,groupby='Age',rotation=90,multi_panel=True,jitter=0)


# #### Visualize in BulkRNA

# In[33]:


visuz.GeneExpression.volcano(df=p10vsp60_combine, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='p10vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.42,1),
                             r=300,show=True,genenames =tuple(P11P60_overlap))


# ## P3 vs P11

# In[34]:


p3vsp11 = sc.get.rank_genes_groups_df(endoth_DE,group='P3',key='wilcoxon_P3vsP11')
p3vsp11 = p3vsp11.set_index(keys='names')
p3vsp11_sig = p3vsp11[(p3vsp11['pvals_adj']<0.05)& (abs(p3vsp11['logfoldchanges'])>2)]


# ### Identify DEs identified in both Bulk and scRNA
# P4 vs P10

# In[35]:


p4vsp10_combine_sig = p4vsp10_combine[(abs(p4vsp10_combine['log2FoldChange'])>2)&(p4vsp10_combine['padj']<0.05)]

P3P11_overlap = []
for gene in p4vsp10_combine_sig['MGI.symbol'].tolist():
    if gene in p3vsp11_sig.index.tolist():
        P3P11_overlap.append(gene)

print('The number of overlapped genes: ',len(P3P11_overlap))


# No overlaped genes in P3 vs P11

# ## Overlap summary
# Order genes by up-regulation and down-regulation in P3 or P11

# In[36]:


#the total overlapped genes
overlap_endoth = P3P60_overlap+P11P60_overlap+P3P11_overlap
overlap_endoth = list(set(overlap_endoth))
print('the total number of DEs identified in both data set:',len(overlap_mesen))

overlap_endoth_df = p3vsp60.loc[overlap_endoth]
overlap_endoth_df=overlap_endoth_df.sort_values('scores',ascending=False)

sc.set_figure_params(scanpy=True, fontsize=14)
sc.pl.heatmap(endoth_DE,var_names=overlap_endoth_df.index.tolist(),groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(12,6))


# ## Unique to scRNA

# In[37]:


p3vsp60_unique = []
for gene in p3vsp60_sig.index.tolist():
    if gene not in overlap_endoth:
        p3vsp60_unique.append(gene)
print('The number of DEs only found in scRNA between P3 vs P60: ',len(p3vsp60_unique))

p11vsp60_unique = []
for gene in p11vsp60_sig.index.tolist():
    if gene not in overlap_endoth:
        p11vsp60_unique.append(gene)
print('The number of DEs only found in scRNA between P11 vs P60: ',len(p11vsp60_unique))

p3vsp11_unique = []
for gene in p3vsp11_sig.index.tolist():
    if gene not in overlap_endoth:
        p3vsp11_unique.append(gene)
print('The number of DEs only found in scRNA between P3 vs P11: ',len(p3vsp11_unique))

#the total unique DEs in scRNA
scrna_DEs_endoth = list(set(p3vsp60_unique+p11vsp60_unique+p3vsp11_unique))
print('The total number of DEs found in scRNA: ',len(scrna_DEs_endoth))
scrna_DEs_endoth_df = p3vsp60.loc[scrna_DEs_endoth]
scrna_DEs_endoth_df=scrna_DEs_endoth_df.sort_values('scores',ascending=False)

sc.set_figure_params(scanpy=True, fontsize=10)
sc.pl.heatmap(endoth_DE,var_names=scrna_DEs_endoth_df.index.tolist(),groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(28,6))


# # Immune

# In[38]:


immune_DE = sc.read_h5ad('../results/Immune_DE.h5ad')


# In[39]:


sc.tl.rank_genes_groups(immune_DE, groupby='Age', groups=['P3','P60'],reference='P60',method = 'wilcoxon', key_added = "wilcoxon_P3vsP60",use_raw=False)
sc.tl.rank_genes_groups(immune_DE, groupby='Age', groups=['P3','P60'],reference='P3',method = 'wilcoxon', key_added = "wilcoxon_P60vsP3",use_raw=False)

sc.tl.rank_genes_groups(immune_DE, groupby='Age', groups=['P11','P60'],reference='P60',method = 'wilcoxon', key_added = "wilcoxon_P11vsP60",use_raw=False)
sc.tl.rank_genes_groups(immune_DE, groupby='Age', groups=['P11','P60'],reference='P11',method = 'wilcoxon', key_added = "wilcoxon_P60vsP11",use_raw=False)

sc.tl.rank_genes_groups(immune_DE, groupby='Age', groups=['P3','P11'],reference='P11',method = 'wilcoxon', key_added = "wilcoxon_P3vsP11",use_raw=False)
sc.tl.rank_genes_groups(immune_DE, groupby='Age', groups=['P3','P11'],reference='P3',method = 'wilcoxon', key_added = "wilcoxon_P11vsP3",use_raw=False)

sc.tl.rank_genes_groups(immune_DE, groupby='Age', method = 'wilcoxon', key_added = "wilcoxon_Age",use_raw=False)


# ## P3 vs P60

# In[40]:


p3vsp60 = sc.get.rank_genes_groups_df(immune_DE,group='P3',key='wilcoxon_P3vsP60')
p3vsp60 =p3vsp60.set_index(keys='names')
p3vsp60_sig = p3vsp60[(p3vsp60['pvals_adj']<0.05)& (abs(p3vsp60['logfoldchanges'])>2)]


# ### DEs identified in both Bulk and scRNA
# P4 vs P60

# In[41]:


p4vsp60_combine_sig = p4vsp60_combine[(abs(p4vsp60_combine['log2FoldChange'])>2)&(p4vsp60_combine['padj']<0.05)]

P3P60_overlap = []
for gene in p4vsp60_combine_sig['MGI.symbol'].tolist():
    if gene in p3vsp60_sig.index.tolist():
        P3P60_overlap.append(gene)

print('The number of overlapped genes: ',len(P3P60_overlap))


# #### Visualize in scRNA 

# In[42]:


for i in np.arange(6,len(P3P60_overlap)+6,6):
    sc.pl.violin(immune_DE,keys=P3P60_overlap[(i-6):i],ncols=6,wspace=0.1,groupby='Age',rotation=90,multi_panel=True,jitter=0)


# #### Visualize in BulkRNA

# In[43]:


visuz.GeneExpression.volcano(df=p4vsp60_combine, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='4vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.42,1),
                             r=300,show=True,genenames =tuple(P3P60_overlap))


# ## P11 vs P60

# In[44]:


p11vsp60 = sc.get.rank_genes_groups_df(immune_DE,group='P11',key='wilcoxon_P11vsP60')
p11vsp60 =p11vsp60.set_index(keys='names')
p11vsp60_sig = p11vsp60[(p11vsp60['pvals_adj']<0.05)& (abs(p11vsp60['logfoldchanges'])>2)]


# ### Identify DEs identified in both Bulk and scRNA
# P10 vs P60

# In[45]:


p10vsp60_combine_sig = p10vsp60_combine[(abs(p10vsp60_combine['log2FoldChange'])>2)&(p10vsp60_combine['padj']<0.05)]

P11P60_overlap = []
for gene in p10vsp60_combine_sig['MGI.symbol'].tolist():
    if gene in p11vsp60_sig.index.tolist():
        P11P60_overlap.append(gene)

print('The number of overlapped genes: ',len(P11P60_overlap))


# #### Visualize in scRNA 

# In[46]:


for i in np.arange(6,len(P11P60_overlap)+6,6):
    sc.pl.violin(immune_DE,keys=P11P60_overlap[(i-6):i],ncols=6,wspace=0.1,groupby='Age',rotation=90,multi_panel=True,jitter=0)


# #### Visualize in BulkRNA

# In[47]:


visuz.GeneExpression.volcano(df=p10vsp60_combine, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='p10vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.42,1),
                             r=300,show=True,genenames =tuple(P11P60_overlap))


# ## P3 vs P11

# In[48]:


p3vsp11 = sc.get.rank_genes_groups_df(immune_DE,group='P3',key='wilcoxon_P3vsP11')
p3vsp11 = p3vsp11.set_index(keys='names')
p3vsp11_sig = p3vsp11[(p3vsp11['pvals_adj']<0.05)& (abs(p3vsp11['logfoldchanges'])>2)]


# ### Identify DEs identified in both Bulk and scRNA
# P4 vs P10

# In[49]:


p4vsp10_combine_sig = p4vsp10_combine[(abs(p4vsp10_combine['log2FoldChange'])>2)&(p4vsp10_combine['padj']<0.05)]

P3P11_overlap = []
for gene in p4vsp10_combine_sig['MGI.symbol'].tolist():
    if gene in p3vsp11_sig.index.tolist():
        P3P11_overlap.append(gene)

print('The number of overlapped genes: ',len(P3P11_overlap))


# No overlaped genes in P3 vs P11

# ## Overlap summary
# Order genes by up-regulation and down-regulation in P3 or P11

# In[50]:


#the total overlapped genes
overlap_immune = P3P60_overlap+P11P60_overlap+P3P11_overlap
overlap_immune = list(set(overlap_immune))
print('the total number of DEs identified in both data set:',len(overlap_immune))

overlap_immune_df = p3vsp60.loc[overlap_immune]
overlap_immune_df=overlap_immune_df.sort_values('scores',ascending=False)
sc.set_figure_params(scanpy=True, fontsize=16)
sc.pl.heatmap(immune_DE,var_names=overlap_immune_df.index.tolist(),groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(12,6))


# ## Unique to scRNA

# In[51]:


p3vsp60_unique = []
for gene in p3vsp60_sig.index.tolist():
    if gene not in overlap_immune:
        p3vsp60_unique.append(gene)
print('The number of DEs only found in scRNA between P3 vs P60: ',len(p3vsp60_unique))

p11vsp60_unique = []
for gene in p11vsp60_sig.index.tolist():
    if gene not in overlap_immune:
        p11vsp60_unique.append(gene)
print('The number of DEs only found in scRNA between P11 vs P60: ',len(p11vsp60_unique))

p3vsp11_unique = []
for gene in p3vsp11_sig.index.tolist():
    if gene not in overlap_immune:
        p3vsp11_unique.append(gene)
print('The number of DEs only found in scRNA between P3 vs P11: ',len(p3vsp11_unique))

#the total unique DEs in scRNA
scrna_DEs_immune = list(set(p3vsp60_unique+p11vsp60_unique+p3vsp11_unique))
print('The total number of DEs found in scRNA: ',len(scrna_DEs_immune))
scrna_DEs_immune_df = p3vsp60.loc[scrna_DEs_immune]
scrna_DEs_immune_df=scrna_DEs_immune_df.sort_values('scores',ascending=False)
sc.set_figure_params(scanpy=True, fontsize=14)
sc.pl.heatmap(immune_DE,var_names=scrna_DEs_immune_df.index.tolist(),groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(14,6))


# # Epithelial

# In[52]:


epith_DE = sc.read_h5ad('../results/Epithelial_DE.h5ad')


# In[53]:


sc.tl.rank_genes_groups(epith_DE, groupby='Age', groups=['P3','P60'],reference='P60',method = 'wilcoxon', key_added = "wilcoxon_P3vsP60",use_raw=False)
sc.tl.rank_genes_groups(epith_DE, groupby='Age', groups=['P3','P60'],reference='P3',method = 'wilcoxon', key_added = "wilcoxon_P60vsP3",use_raw=False)

sc.tl.rank_genes_groups(epith_DE, groupby='Age', groups=['P11','P60'],reference='P60',method = 'wilcoxon', key_added = "wilcoxon_P11vsP60",use_raw=False)
sc.tl.rank_genes_groups(epith_DE, groupby='Age', groups=['P11','P60'],reference='P11',method = 'wilcoxon', key_added = "wilcoxon_P60vsP11",use_raw=False)

sc.tl.rank_genes_groups(epith_DE, groupby='Age', groups=['P3','P11'],reference='P11',method = 'wilcoxon', key_added = "wilcoxon_P3vsP11",use_raw=False)
sc.tl.rank_genes_groups(epith_DE, groupby='Age', groups=['P3','P11'],reference='P3',method = 'wilcoxon', key_added = "wilcoxon_P11vsP3",use_raw=False)

sc.tl.rank_genes_groups(epith_DE, groupby='Age', method = 'wilcoxon', key_added = "wilcoxon_Age",use_raw=False)


# ## P3 vs P60

# In[54]:


p3vsp60 = sc.get.rank_genes_groups_df(epith_DE,group='P3',key='wilcoxon_P3vsP60')
p3vsp60 =p3vsp60.set_index(keys='names')
p3vsp60_sig = p3vsp60[(p3vsp60['pvals_adj']<0.05)& (abs(p3vsp60['logfoldchanges'])>2)]


# ### DEs identified in both Bulk and scRNA
# P4 vs P60

# In[55]:


p4vsp60_combine_sig = p4vsp60_combine[(abs(p4vsp60_combine['log2FoldChange'])>2)&(p4vsp60_combine['padj']<0.05)]

P3P60_overlap = []
for gene in p4vsp60_combine_sig['MGI.symbol'].tolist():
    if gene in p3vsp60_sig.index.tolist():
        P3P60_overlap.append(gene)

print('The number of overlapped genes: ',len(P3P60_overlap))


# #### Visualize in scRNA 

# In[56]:


for i in np.arange(6,len(P3P60_overlap)+6,6):
    sc.pl.violin(epith_DE,keys=P3P60_overlap[(i-6):i],ncols=6,wspace=0.1,groupby='Age',rotation=90,multi_panel=True,jitter=0)


# #### Visualize in BulkRNA

# In[57]:


visuz.GeneExpression.volcano(df=p4vsp60_combine, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='4vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.42,1),
                             r=300,show=True,genenames =tuple(P3P60_overlap[:40]))


# ## P11 vs P60

# In[58]:


p11vsp60 = sc.get.rank_genes_groups_df(epith_DE,group='P11',key='wilcoxon_P11vsP60')
p11vsp60 =p11vsp60.set_index(keys='names')
p11vsp60_sig = p11vsp60[(p11vsp60['pvals_adj']<0.05)& (abs(p11vsp60['logfoldchanges'])>2)]


# ### Identify DEs identified in both Bulk and scRNA
# P10 vs P60

# In[59]:


p10vsp60_combine_sig = p10vsp60_combine[(abs(p10vsp60_combine['log2FoldChange'])>2)&(p10vsp60_combine['padj']<0.05)]

P11P60_overlap = []
for gene in p10vsp60_combine_sig['MGI.symbol'].tolist():
    if gene in p11vsp60_sig.index.tolist():
        P11P60_overlap.append(gene)

print('The number of overlapped genes: ',len(P11P60_overlap))


# #### Visualize in scRNA 

# In[60]:


for i in np.arange(6,len(P11P60_overlap)+6,6):
    sc.pl.violin(epith_DE,keys=P11P60_overlap[(i-6):i],ncols=6,wspace=0.1,groupby='Age',rotation=90,multi_panel=True,jitter=0)


# #### Visualize in BulkRNA

# In[61]:


visuz.GeneExpression.volcano(df=p10vsp60_combine, lfc='log2FoldChange', pv='pvalue', dotsize=0.5,geneid='GeneNames',sign_line=True, axtickfontsize=10,
                             axtickfontname='Verdana',gfont=6,figname='p10vsp60',dim=(6,6),plotlegend=True,legendpos='upper right',legendanchor=(1.42,1),
                             r=300,show=True,genenames =tuple(P11P60_overlap))


# ## P3 vs P11

# In[62]:


p3vsp11 = sc.get.rank_genes_groups_df(epith_DE,group='P3',key='wilcoxon_P3vsP11')
p3vsp11 = p3vsp11.set_index(keys='names')
p3vsp11_sig = p3vsp11[(p3vsp11['pvals_adj']<0.05)& (abs(p3vsp11['logfoldchanges'])>2)]


# ### Identify DEs identified in both Bulk and scRNA
# P4 vs P10

# In[63]:


p4vsp10_combine_sig = p4vsp10_combine[(abs(p4vsp10_combine['log2FoldChange'])>2)&(p4vsp10_combine['padj']<0.05)]

P3P11_overlap = []
for gene in p4vsp10_combine_sig['MGI.symbol'].tolist():
    if gene in p3vsp11_sig.index.tolist():
        P3P11_overlap.append(gene)

print('The number of overlapped genes: ',len(P3P11_overlap))


# No overlaped genes in P3 vs P11

# ## Overlap summary
# Order genes by up-regulation and down-regulation in P3 or P11

# In[64]:


#the total overlapped genes
overlap_epith = P3P60_overlap+P11P60_overlap+P3P11_overlap
overlap_epith = list(set(overlap_epith))
print('the total number of DEs identified in both data set:',len(overlap_epith))

overlap_epith_df = p3vsp60.loc[overlap_epith]
overlap_epith_df=overlap_epith_df.sort_values('scores',ascending=False)
sc.set_figure_params(scanpy=True, fontsize=10)
sc.pl.heatmap(epith_DE,var_names=overlap_epith_df.index.tolist(),groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(24,6))


# ## Unique to scRNA

# In[65]:


p3vsp60_unique = []
for gene in p3vsp60_sig.index.tolist():
    if gene not in overlap_epith:
        p3vsp60_unique.append(gene)
print('The number of DEs only found in scRNA between P3 vs P60: ',len(p3vsp60_unique))

p11vsp60_unique = []
for gene in p11vsp60_sig.index.tolist():
    if gene not in overlap_epith:
        p11vsp60_unique.append(gene)
print('The number of DEs only found in scRNA between P11 vs P60: ',len(p11vsp60_unique))

p3vsp11_unique = []
for gene in p3vsp11_sig.index.tolist():
    if gene not in overlap_epith:
        p3vsp11_unique.append(gene)
print('The number of DEs only found in scRNA between P3 vs P11: ',len(p3vsp11_unique))

#the total unique DEs in scRNA
scrna_DEs_epith = list(set(p3vsp60_unique+p11vsp60_unique+p3vsp11_unique))
print('The total number of DEs found in scRNA: ',len(scrna_DEs_epith))
scrna_DEs_epith_df = p3vsp60.loc[scrna_DEs_epith]
scrna_DEs_epith_df=scrna_DEs_epith_df.sort_values('scores',ascending=False)

sc.pl.heatmap(epith_DE,var_names=scrna_DEs_epith_df.index.tolist(),groupby="Age", show_gene_labels=False,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(20,6))


# ### Up-regulated genes

# In[66]:


up_gene = scrna_DEs_epith_df[scrna_DEs_epith_df.logfoldchanges>0].index.tolist()
sc.set_figure_params(scanpy=True, fontsize=10)
sc.pl.heatmap(epith_DE,var_names=up_gene[:180],groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(26,6))
sc.pl.heatmap(epith_DE,var_names=up_gene[180:],groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(26,6))


# ### Down-regulated genes

# In[67]:


down_gene_df = scrna_DEs_epith_df[scrna_DEs_epith_df.logfoldchanges<0]
down_gene_df = down_gene_df.sort_values('scores',ascending=True)
down_gene = down_gene_df.index.to_list()

sc.set_figure_params(scanpy=True, fontsize=10)
sc.pl.heatmap(epith_DE,var_names=down_gene[:160],groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(26,6))
sc.pl.heatmap(epith_DE,var_names=down_gene[160:],groupby="Age", show_gene_labels=True,use_raw=False,
                                layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(26,6))


# Highlight each 50th gene
# ax_dict = sc.pl.heatmap(epith_DE,var_names=scrna_DEs_epith_df.index.tolist(),groupby="Age", show_gene_labels=True,use_raw=False,show=False,
#                                 layer='scaled',dendrogram=False,vmax=-5, vmin=5, cmap='bwr',figsize=(20,6))
#                                 
# ax_dict['heatmap_ax'].set_xticks(range(len(scrna_DEs_epith_df.index.tolist()))[::50])
# 
# ax_dict['heatmap_ax'].set_xticklabels(scrna_DEs_epith_df.index.tolist()[::50])               
# 

# # Save the unique and overlapped DEs 
# 
# subset unique/overlapped genes in P3 vs P60 in all four cell types

# In[174]:


overlap_mesen_df.to_csv('OverlappedDEs_Mesenchymal_BulkRNA.csv')
scrna_DEs_mesen_df.to_csv('UniqueDEs_Mesenchymal.csv')

overlap_endoth_df.to_csv('OverlappedDEs_Endothelial_BulkRNA.csv')
scrna_DEs_endoth_df.to_csv('UniqueDEs_Endothelial.csv')

overlap_immune_df.to_csv('OverlappedDEs_Immune_BulkRNA.csv')
scrna_DEs_immune_df.to_csv('UniqueDEs_Immune.csv')

overlap_epith_df.to_csv('OverlappedDEs_Epithelial_BulkRNA.csv')
scrna_DEs_epith_df.to_csv('UniqueDEs_Epithelial.csv')


# # Compare overlapped/unique DEs across all four cell types

# In[13]:


adata = sc.read_h5ad('../../results/Mouse_harmony_rmBatchbetweenSamples_rmDoublets.h5ad')


# In[14]:


adata1 = adata[(adata.obs['CellType_edited4']!='NeuroGlia')& ((adata.obs['Age']!='Adult'))]
adata1 = adata1.raw.to_adata()
sc.pp.scale(adata1)


# Define heatmap function

# In[12]:


from itertools import groupby
import datetime

def add_line(ax, xpos, ypos):
    line = plt.Line2D([ypos, ypos+ .2], [xpos, xpos], color='black', transform=ax.transAxes)
    line.set_clip_on(False)
    ax.add_line(line)

def label_len(my_index,level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k,g in groupby(labels)]

def label_group_bar_table(ax, df):
    xpos = -.2
    scale = 1./df.index.size
    for level in range(df.index.nlevels):
        pos = df.index.size
        for label, rpos in label_len(df.index,level):
            add_line(ax, pos*scale, xpos)
            pos -= rpos
            lypos = (pos + .5 * rpos)*scale
            ax.text(xpos+.1, lypos, label, ha='center', transform=ax.transAxes) 
        add_line(ax, pos*scale , xpos)
        xpos -= .2


# ## Mesenchymal

# ### overlapped DEs only in Mesenchymal

# In[319]:


overlap_gene_list1 = []
for gene in overlap_mesen:
    if (gene not in overlap_endoth) and (gene not in overlap_immune) and (gene not in overlap_epith):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)


# In[320]:


subset_adata = adata1[:, overlap_gene_list1].copy()
subset = adata1[:, overlap_gene_list1].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[321]:


fig = plt.figure(figsize = (10, 14))
ax = fig.add_subplot(111)
sns.heatmap(subset,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset)
fig.subplots_adjust(bottom=.1*subset.index.nlevels)
plt.show()


# ### Overlapped DEs shared among cell types

# In[298]:


overlap_gene_list1 = []
for gene in overlap_mesen:
    if (gene in overlap_endoth) or (gene in overlap_immune) or (gene in overlap_epith):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)


# In[299]:


subset_adata = adata1[:, overlap_gene_list1].copy()
subset = adata1[:, overlap_gene_list1].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[300]:


fig = plt.figure(figsize = (10, 14))
ax = fig.add_subplot(111)
sns.heatmap(subset,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset)
fig.subplots_adjust(bottom=.1*subset.index.nlevels)
plt.show()


# ## Endothelial

# ### Overlapped DEs only in Endothelial

# In[322]:


overlap_gene_list1 = []
for gene in overlap_endoth:
    if (gene not in overlap_mesen) and (gene not in overlap_immune) and (gene not in overlap_epith):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)


# In[323]:


subset_adata = adata1[:, overlap_gene_list1].copy()
subset = adata1[:, overlap_gene_list1].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[324]:


fig = plt.figure(figsize = (12, 14))
ax = fig.add_subplot(111)
sns.heatmap(subset,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset)
fig.subplots_adjust(bottom=.1*subset.index.nlevels)
plt.show()


# ### Overlapped DEs shared among cell types

# In[325]:


overlap_gene_list1 = []
for gene in overlap_endoth:
    if (gene in overlap_mesen) or (gene in overlap_immune) or (gene in overlap_epith):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)


# In[305]:


subset_adata = adata1[:, overlap_gene_list1].copy()
subset = adata1[:, overlap_gene_list1].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[307]:


fig = plt.figure(figsize = (8, 14))
ax = fig.add_subplot(111)
sns.heatmap(subset,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset)
fig.subplots_adjust(bottom=.1*subset.index.nlevels)
plt.show()


# ## Immune

# ### overlapped DEs only in Immune cells

# In[286]:


overlap_gene_list1 = []
for gene in overlap_immune:
    if (gene not in overlap_mesen) and (gene not in overlap_endoth) or (gene not in overlap_epith):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)


# In[287]:


subset_adata = adata1[:, overlap_gene_list1].copy()
subset = adata1[:, overlap_gene_list1].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[288]:


fig = plt.figure(figsize = (8, 16))
ax = fig.add_subplot(111)
sns.heatmap(subset,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset)
fig.subplots_adjust(bottom=.1*subset.index.nlevels)
plt.show()


# ### Overlapped DEs shared among cell types

# In[316]:


overlap_gene_list1 = []
for gene in overlap_immune:
    if (gene in overlap_mesen) or (gene in overlap_endoth) or (gene in overlap_epith):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)


# In[317]:


subset_adata = adata1[:, overlap_gene_list1].copy()
subset = adata1[:, overlap_gene_list1].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[318]:


fig = plt.figure(figsize = (8, 16))
ax = fig.add_subplot(111)
sns.heatmap(subset,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset)
fig.subplots_adjust(bottom=.1*subset.index.nlevels)
plt.show()


# ## Epithelial

# ### Overlapped DEs only in Epithelial 

# In[289]:


overlap_gene_list1 = []
for gene in overlap_epith:
    if (gene not in overlap_mesen) and (gene not in overlap_endoth) and (gene not in overlap_immune):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)


# In[290]:


subset_adata = adata1[:, overlap_gene_list1].copy()
subset = adata1[:, overlap_gene_list1].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[291]:


len(overlap_gene_list1)


# In[292]:


subset1 = subset.iloc[:,:73]
fig = plt.figure(figsize = (24, 16))
ax = fig.add_subplot(111)
sns.heatmap(subset1,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset1)
fig.subplots_adjust(bottom=.1*subset1.index.nlevels)
plt.show()


# In[293]:


subset2 = subset.iloc[:,73:len(overlap_gene_list1)]
fig = plt.figure(figsize = (24, 16))
ax = fig.add_subplot(111)
sns.heatmap(subset2,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset2)
fig.subplots_adjust(bottom=.1*subset2.index.nlevels)
plt.show()


# ### Overlapped DEs shared among cells types

# In[312]:


overlap_gene_list1 = []
for gene in overlap_epith:
    if (gene in overlap_mesen) and (gene in overlap_endoth) or (gene in overlap_immune):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)


# In[313]:


subset_adata = adata1[:, overlap_gene_list1].copy()
subset = adata1[:, overlap_gene_list1].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[314]:


len(overlap_gene_list1)


# In[315]:


fig = plt.figure(figsize = (8, 16))
ax = fig.add_subplot(111)
sns.heatmap(subset,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset)
fig.subplots_adjust(bottom=.1*subset.index.nlevels)
plt.show()


# # Visualize all overlapped DEs among cell types

# In[15]:


overlap_mesen_df = pd.read_csv('OverlappedDEs_Mesenchymal_BulkRNA.csv')

overlap_endoth_df = pd.read_csv('OverlappedDEs_Endothelial_BulkRNA.csv')

overlap_immune_df = pd.read_csv('OverlappedDEs_Immune_BulkRNA.csv')

overlap_epith_df = pd.read_csv('OverlappedDEs_Epithelial_BulkRNA.csv')


# In[16]:


overlap_mesen = overlap_mesen_df.names.tolist()
overlap_endoth = overlap_endoth_df.names.tolist()
overlap_immune = overlap_immune_df.names.tolist()
overlap_epith = overlap_epith_df.names.tolist()


# In[17]:


#mesenchymal
overlap_gene_list1 = []
for gene in overlap_mesen:
    if (gene in overlap_endoth) or (gene in overlap_immune) or (gene in overlap_epith):
        overlap_gene_list1.append(gene)
print(overlap_gene_list1)

#Endothelial
overlap_gene_list2 = []
for gene in overlap_endoth:
    if (gene in overlap_mesen) or (gene in overlap_immune) or (gene in overlap_epith):
        overlap_gene_list2.append(gene)
print(overlap_gene_list2)

#Immune
overlap_gene_list3 = []
for gene in overlap_immune:
    if (gene not in overlap_mesen) or (gene not in overlap_endoth) or (gene  in overlap_epith):
        overlap_gene_list3.append(gene)
print(overlap_gene_list3)

#epithelial
overlap_gene_list4 = []
for gene in overlap_epith:
    if (gene in overlap_mesen) or (gene in overlap_endoth) or (gene in overlap_immune):
        overlap_gene_list4.append(gene)
print(overlap_gene_list4)


# In[28]:


overlap_gene_list = overlap_gene_list1+overlap_gene_list2+overlap_gene_list3+overlap_gene_list4
used = []
unique = [used.append(x) for x in overlap_gene_list if x not in used]
print (len(unique))


# In[30]:


subset_adata = adata1[:, used].copy()
subset = adata1[:, used].to_df()
subset= subset.set_index([subset_adata.obs['CellType_edited4'],subset_adata.obs['Age']])
subset = subset.sort_index(level=0)


# In[32]:


fig = plt.figure(figsize = (20, 14))
ax = fig.add_subplot(111)
sns.heatmap(subset,vmax=-5, vmin=5, cmap='bwr')
#Below 3 lines remove default labels
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')

label_group_bar_table(ax, subset)
fig.subplots_adjust(bottom=.1*subset.index.nlevels)
plt.show()


# In[39]:


data = {'Genes': used}
df = pd.DataFrame(data) 
df.to_csv('DEs_Bulk_scRNA_shared_allCellTypes.csv')

