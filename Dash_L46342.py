
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import skimage
import tifffile
import skimage.feature as ft
import scipy.ndimage as ndimage

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input,Output
import plotly.express as px
import scanpy as sc

data='L46342v2/'

hvadata = sc.read(data+'hvadata.h5ad')
adata = sc.read(data+'adata.h5ad')

pos_df=pd.read_csv(data+'masked_coords.csv')
pos_df['IntX']=np.floor(pos_df['x']*1).astype(int)
pos_df['IntY']=np.floor(pos_df['y']*1).astype(int)
#pos_df['Counts']=np.sum(gene_img, axis=1)
pos_df['Counts']=adata.obs['total_counts'].values

adata.xpos = adata.obs.merge(pos_df, left_index=True, right_on = 'cells')['x'].values
adata.ypos = adata.obs.merge(pos_df, left_index=True, right_on = 'cells')['y'].values

hvadata.xpos = hvadata.obs.merge(pos_df, left_index=True, right_on = 'cells')['x'].values
hvadata.ypos = hvadata.obs.merge(pos_df, left_index=True, right_on = 'cells')['y'].values

df = pd.DataFrame({'UMAPx':hvadata.obsm['X_umap'][:,0], 'UMAPy':hvadata.obsm['X_umap'][:,1], 
                   'X':hvadata.xpos, 'Y':hvadata.ypos, 'Cluster':hvadata.obs['leiden'].values,
                  'Counts':hvadata.obs['total_counts'].values, 'Unnamed: 0':hvadata.obs.index})

real_vals = sp.sparse.csr_matrix.toarray(adata.X)
smoothed_data_02 = np.load(data+'smoothed_data_02.npy')

def get_base_UMAP_trace():
    traces=[]
    for c in df['Cluster'].unique():
        ddf = df[df['Cluster']==c]
        traces.append(go.Scattergl(x=ddf['UMAPx'],y=ddf['UMAPy'],mode='markers',text=ddf.index.values, 
                               name='Cluster '+c, marker=dict(size=4))) #visible='legendonly'
    layout=go.Layout(title='UMAP', dragmode='lasso')
    return {'data':traces,'layout':layout}

def get_filtered_UMAP_trace(cluster):
    traces=[]
    for c in df['Cluster'].unique():
        ddf = df[df['Cluster']==c]
        if c==cluster:
            traces.append(go.Scattergl(x=ddf['UMAPx'],y=ddf['UMAPy'],mode='markers',text=ddf.index.values, 
                               name='Cluster '+c, marker=dict(size=4))) 
        else:
            traces.append(go.Scattergl(x=ddf['UMAPx'],y=ddf['UMAPy'],mode='markers',text=ddf.index.values, 
                               name='Cluster '+c, visible='legendonly', marker=dict(size=4))) 
    layout=go.Layout(title='UMAP', dragmode='lasso')
    return {'data':traces,'layout':layout}

def get_base_slideseq_trace():
    trace1 = go.Scattergl(x=df['Y'],y=-df['X'],marker=dict(size=df['Counts']/500), mode='markers')
    layout=go.Layout(title="Bead Hits",height=good_height,template='plotly_dark')
    return {'data':[trace1], 'layout':layout} 

def get_cluster_filtered_slideseq_trace(cluster):
    f_df=df[df['Cluster']==cluster]
    return get_filtered_slideseq_trace(f_df)

def get_index_filtered_slideseq_trace(indices):
    f_df=df.loc[indices]
    return get_filtered_slideseq_trace(f_df)

def get_filtered_slideseq_trace(f_df):
    trace2 = go.Scattergl(x=df['Y'],y=-df['X'],marker=dict(size=df['Counts']/500), mode='markers')
    trace3 = go.Scattergl(x=f_df['Y'],y=-f_df['X'],marker=dict(size=8), mode='markers')
    layout2=layout=go.Layout(title="Bead Hits",height=good_height,template='plotly_dark')

    return {'data':[trace2, trace3],'layout':layout2}  

def get_raw_gene_image(gene):
    idx = np.where(adata.var.index==gene)[0][0]
    tmp=real_vals[:,idx]
    sizes=2*(tmp-np.min(tmp))/(np.max(tmp)-np.min(tmp))
    trace2 = go.Scattergl(x=df['Y'],y=-df['X'],marker=dict(size=df['Counts']/500), mode='markers', hoverinfo='none')
    trace3 = go.Scattergl(x=df['Y'],y=-df['X'],marker=dict(size=15*sizes), mode='markers',
                             hovertext='')
    layout2=layout=go.Layout(title="Raw UMIs:  "+gene,height=good_height,template='plotly_dark')
    return {'data':[trace2, trace3],'layout':layout2}  

def get_gene_image(gene):
    idx = np.where(adata.var.index==gene)[0][0]
    tmp=smoothed_data_02[:,idx]
    sizes=2*(tmp-np.min(tmp))/(np.max(tmp)-np.min(tmp))
    trace2 = go.Scattergl(x=df['Y'],y=-df['X'],marker=dict(size=df['Counts']/500), mode='markers', hoverinfo='none')
    trace3 = go.Scattergl(x=df['Y'],y=-df['X'],marker=dict(size=15*sizes), mode='markers',
                             hovertext='')
    layout2=layout=go.Layout(title="Smoothed UMIs:  "+gene,height=good_height,template='plotly_dark')
    return {'data':[trace2, trace3],'layout':layout2}  

def get_gene_list(indices):
    log2ratio=np.mean(smoothed_data_02[indices], axis=0)/np.mean(smoothed_data_02, axis=0)
    log2ratio=log2ratio*(np.mean(smoothed_data_02[indices], axis=0)>0.35).astype(int)
    log2idx=np.argsort(-log2ratio)
        
    return [{'label': adata.var.index[log2idx[i]], 'value': adata.var.index[log2idx[i]]} for i in range(0,20)]

base_gene_list = [
            {'label': 'ESR2', 'value': 'ESR2'},
            
        ]
default_gene = 'ESR2'



# Build simple Dash Layout

app = dash.Dash(__name__)

good_height=800

cluster_list = [{'label':'All', 'value':'All'}]
cluster_list = cluster_list+[{'label':'Cluster '+c, 'value':'Cluster '+c} for c in df['Cluster'].unique()]

app.layout = html.Div([
    dcc.Dropdown(
        id='cluster-dropdown',
        options=cluster_list,
        value='All'
    ),
    dcc.Graph(id='UMAP', figure={'layout': {'height':good_height}}, style={'width': '49%', 'display': 'inline-block'}), #Graph that displays all data
    dcc.Graph(id='ImageGraph', figure={'layout': {'height':good_height}}, style={'width': '49%', 'display': 'inline-block'}), #Graph that shows only filtered data    
    html.Div(id='text_display'),  #To show format of selectData
    html.Button(id='ignore',style={'display':'none'}), #Create a hidden button just so the first callback will have an input.  It doesn't so anything.
    dcc.Dropdown(
        id='gene-dropdown',
        options=base_gene_list,
        value=default_gene
    ),
    dcc.Graph(id='GeneGraph', figure={'layout': {'height':good_height}}, style={'width': '49%', 'display': 'inline-block'}), #Graph that shows only filtered data    
    dcc.Graph(id='RawGeneGraph', figure={'layout': {'height':good_height}}, style={'width': '49%', 'display': 'inline-block'}), #Graph that shows only filtered data    
    dcc.Input(id='GeneText', placeholder='Enter a gene...', type='text', value='')
    
])

current_list=[]

#This sets up the callback for the cluster selection dropdown and has it update the UMAP
@app.callback(Output('UMAP', 'figure'), Input('cluster-dropdown', 'value'))
def update_cluster(text_value):
    ctx=dash.callback_context
    if (not ctx.triggered) | (text_value in 'All'):
        return get_base_UMAP_trace()
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    cluster = text_value.split(' ')[1]
    return get_filtered_UMAP_trace(cluster)
    

#This sets up the callback for drawing on the UMAP graph and updating the Slideseq image
@app.callback(Output('ImageGraph','figure'),[Input('UMAP','selectedData'), Input('cluster-dropdown', 'value')])
def selectData3(selectData, text_value):
    
    ctx=dash.callback_context
    if (not ctx.triggered):
        return get_base_slideseq_trace()
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id=='UMAP':   
        if (selectData!=None):
            indices=np.asarray([ f['text'] for f in selectData['points']])

            return get_index_filtered_slideseq_trace(indices)
        else:
            return get_base_slideseq_trace()
    else:
        if text_value in 'All':
            return get_base_slideseq_trace()
        else:
            cluster = text_value.split(' ')[1]

            return get_index_filtered_slideseq_trace(df['Cluster']==cluster)
        

@app.callback(Output('gene-dropdown', 'options'), [Input('UMAP', 'selectedData'), Input('cluster-dropdown', 'value')])
def update_date_dropdown(selectData, text_value):
    ctx=dash.callback_context
    if (not ctx.triggered):
        return base_gene_list
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id=='UMAP':   
        if (selectData!=None):    
            indices=np.asarray([ f['text'] for f in selectData['points']])
            return get_gene_list(indices)

        else:
            return base_gene_list        
    else:
        if text_value in 'All':
            return base_gene_list
        else:
            cluster = text_value.split(' ')[1]
            return get_gene_list(df['Cluster']==cluster)

@app.callback(Output('RawGeneGraph','figure'),[Input('gene-dropdown', 'value'), Input('GeneText', 'value')])
def selectData3(selectGene, textvalue):    
    ctx=dash.callback_context
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if not ctx.triggered:
        return get_raw_gene_image(default_gene)
    if trigger_id=='GeneText':
        if (len(np.where(adata.var.index==textvalue)[0])>0):
            return get_raw_gene_image(textvalue)
        else:
            return get_raw_gene_image(default_gene)
    else:
        return get_raw_gene_image(selectGene)
    


@app.callback(Output('GeneGraph','figure'),[Input('gene-dropdown', 'value'), Input('GeneText', 'value')])
def selectData4(selectGene, textvalue):    
    ctx=dash.callback_context
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if not ctx.triggered:
        return get_gene_image(default_gene)
    
    if trigger_id=='GeneText':
        if (len(np.where(adata.var.index==textvalue)[0])>0):
            return get_gene_image(textvalue)
        else:
            return get_gene_image(default_gene)
    else:
        return get_gene_image(selectGene)

#Dropdown response:  Raw Gene Graph







if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8092, debug=True) #To let others connect
    #app.run_server(debug=True, use_reloader=False) #For debugging