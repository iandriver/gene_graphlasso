#!/usr/bin/python3
import sys
import datetime
import numpy as np
from sklearn.covariance import GraphLasso
from sklearn.decomposition import PCA as skPCA
import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg')
import itertools
from sklearn.preprocessing import scale
import community
import matplotlib.pyplot as plt
import os
import networkx as nx
from sklearn import cluster, covariance, manifold
from matplotlib import collections
import seaborn as sns

def return_top_pca_gene(by_cell_matrix, range_genes = None):
    gene_number = 100
    gene_pca = skPCA(n_components=3)
    np_by_gene = np.asarray(by_cell_matrix.transpose())
    gene_index = by_cell_matrix.index.tolist()

    if range_genes is not None:
        start_num= range_genes[0]
        end_num_genes = range_genes[1]
    else:
        start_num = 0
        end_num_genes = min(gene_number, len(gene_index))
    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=gene_index)
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len(gene_index))
    top_pca_list = Pc_sort_df.index.tolist()
    new_cell_matrix = by_cell_matrix.ix[top_pca_list[start_num:end_num_genes],:]
    return new_cell_matrix.transpose(), top_pca_list[start_num:end_num_genes]

def save_network_graph(matrix, labels, filename, title, scale=8, node_weight = None, layout = "circular", weight = lambda x: abs(4*x)**(2.5)):
    labels = dict( zip( range( len(labels) ), labels) )
    d = matrix.shape[0]
    D = nx.Graph(matrix)
    weights = [ D[x][y]['weight'] for x,y in D.edges() ]
    large_cutoff = np.mean(weights) + np.std(weights)
    small_cutoff = np.mean(weights) - np.std(weights)
    elarge=[(u,v) for (u,v,d) in D.edges(data=True) if d['weight'] >large_cutoff]
    esmall=[(u,v) for (u,v,d) in D.edges(data=True) if d['weight'] <=small_cutoff]
    #weights = weights/np.max( np.abs( weights ) )
    cmap = plt.get_cmap( "Reds" ) #or some other one

    fig = plt.figure(figsize=(50,50))
    ax = fig.add_subplot(111)
    if layout == "circular":
    	pos = nx.circular_layout( D, scale =scale )
    elif layout == "spring":
    	pos = nx.spring_layout( D ,scale = scale, iterations = 35 )

    bweights = [ 'k'*(z<0) + 'r'*(z>0) for z in weights ]
    width_small = [ weight(w) for w in weights if w <= small_cutoff]
    width_large = [ weight(w) for w in weights if w > large_cutoff]
    nx.draw_networkx_edges(D,pos,edgelist=elarge,
                    edge_color= "green", width=width_large, )
    nx.draw_networkx_edges(D,pos,edgelist=esmall,width=width_small,alpha=0.5,edge_color="red",style='dashed')
    if node_weight == None:
        nx.draw_networkx_nodes( D, pos, ax=ax, node_size = 1, node_color="black")
    else:
        nx.draw_networkx_nodes( D, pos, ax=ax, node_size = node_weight, node_color="black", alpha=0.4)
    nx.draw_networkx_labels( D, pos,font_size=18, labels = labels, ax = ax)
    plt.axis("off")
    plt.title(title)
    plt.savefig(filename, bbox_inches="tight")

def make_pcor_network(qt5_data):
    path_to_file = qt5_data.filepath.text()

    #load file gene
    by_cell = pd.DataFrame.from_csv(path_to_file, sep='\t')
    by_gene = by_cell.transpose()
    #create list of genes
    gene_list = by_cell.index.tolist()
    #create cell list
    cell_list = [x for x in list(by_cell.columns.values)]


    df_by_gene1 = pd.DataFrame(by_gene, columns=gene_list, index=cell_list)
    df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list)
    log2_df_cell = np.log2(df_by_cell1+1)
    alpha=0.4
    keep_genes = []
    num_iterations = int(qt5_data.iterations.text())
    gene_index_list = list(range(0,((num_iterations+1)*100), 100))

    for iter_genes in range(0,num_iterations):
        top_pca_matrix, top_pca_genes = return_top_pca_gene(log2_df_cell,range_genes=[gene_index_list[max(0,iter_genes-1)],gene_index_list[iter_genes+1]])
        mean_expr_dict = {}

        gl = covariance.GraphLassoCV()
        gene_data = scale(top_pca_matrix.as_matrix())

        gl.fit(gene_data)
        _, labels = cluster.affinity_propagation(gl.covariance_)

        n_labels = labels.max()
        names = np.array(top_pca_genes)
        prec_sp = gl.precision_
        matrix1 = -prec_sp + np.diag( np.diagonal( prec_sp) )
        D = nx.Graph(matrix1)

        gene_weights_dict ={}
        for n in names:
            gene_weights_dict[n] = 0
        for x,y in D.edges():
            gene1 = names[x]
            gene2 = names[y]
            abs_weight = abs(D[x][y]['weight'])
            gene_weights_dict[gene1] += abs_weight
            gene_weights_dict[gene2] += abs_weight

        clust_gene_list = []
        avg_wieght_list = []
        for i in range(n_labels + 1):
            clust_id = "cluster_"+str(i+1)
            w_list = [gene_weights_dict[g] for g in names[labels == i]]
            clust_gene_list.append([names[labels == i]])
            avg_wieght_list.append(sum(w_list)/len(w_list))
            print('Cluster %i: %s' % ((i + 1), ', '.join(names[labels == i])))
        if iter_genes == 0:
            threshold = np.mean(avg_wieght_list)-np.std(avg_wieght_list)

        for g_list in [clust_gene_list[i] for i,av in enumerate(avg_wieght_list) if av >= threshold]:
            keep_genes.append(np.ravel(g_list))


    final_gene_list = list(set([item for sublist in keep_genes for item in sublist]))
    print(final_gene_list)
    gene_matrix_log2 = log2_df_cell.T
    top_matrix = gene_matrix_log2[final_gene_list]
    top_matrix_cell = top_matrix.T
    return top_matrix, top_matrix_cell, final_gene_list

def find_node_weight(top_matrix, final_gene_list):
    mean_expr_dict = {}
    for g in final_gene_list:
        g_expr = top_matrix[g]
        mean_expr_dict[g] = np.mean(g_expr)

    node_weights =[int(500*round(v)) for k,v in mean_expr_dict.items()]
    return node_weights

def make_graphlasso(top_matrix, final_gene_list):
    gl = covariance.GraphLassoCV()
    gene_data = scale(top_matrix.as_matrix())

    gl.fit(gene_data)

    names = np.array(final_gene_list)

    node_position_model = manifold.LocallyLinearEmbedding(
        n_components=2, eigen_solver='dense', n_neighbors=7)

    embedding = node_position_model.fit_transform(gene_data.T).T

    return gl, embedding, names

def plot_networks(gl, embedding, names, node_weights, qt5_data):
    _, labels = cluster.affinity_propagation(gl.covariance_)
    n_labels = labels.max()
    path_to_file = qt5_data.filepath.text()
    plt.figure(1, facecolor='w', figsize=(10, 8))
    plt.clf()
    ax = plt.axes([0., 0., 1., 1.])
    plt.axis('off')

    # Display a graph of the partial correlations
    partial_correlations = gl.precision_.copy()
    d = 1 / np.sqrt(np.diag(partial_correlations))
    partial_correlations *= d
    partial_correlations *= d[:, np.newaxis]
    non_zero = (np.abs(np.triu(partial_correlations, k=1)) > 0.02)

    # Plot the nodes using the coordinates of our embedding
    cm = plt.cm.get_cmap('spectral')
    #cmap = [cm(1.*i/len(embedding[0])) for i in range(len(embedding[0]))]


    plt.scatter(embedding[0], embedding[1], s=100 * d ** 2, c=labels,
                cmap=cm)

    # Plot the edges
    start_idx, end_idx = np.where(non_zero)
    #a sequence of (*line0*, *line1*, *line2*), where::
    #            linen = (x0, y0), (x1, y1), ... (xm, ym)
    segments = [[embedding[:, start], embedding[:, stop]]
                for start, stop in zip(start_idx, end_idx)]
    values = np.abs(partial_correlations[non_zero])
    lc = collections.LineCollection(segments,
                        zorder=0, cmap=plt.cm.hot_r,
                        norm=plt.Normalize(0, .7 * values.max()))
    lc.set_array(values)
    lc.set_linewidths(15 * values)
    ax.add_collection(lc)

    # Add a label to each node. The challenge here is that we want to
    # position the labels to avoid overlap with other labels
    for index, (name, label, (x, y)) in enumerate(
            zip(names, labels, embedding.T)):

        dx = x - embedding[0]
        dx[index] = 1
        dy = y - embedding[1]
        dy[index] = 1
        this_dx = dx[np.argmin(np.abs(dy))]
        this_dy = dy[np.argmin(np.abs(dx))]
        if this_dx > 0:
            horizontalalignment = 'left'
            x = x + .002
        else:
            horizontalalignment = 'right'
            x = x - .002
        if this_dy > 0:
            verticalalignment = 'bottom'
            y = y + .002
        else:
            verticalalignment = 'top'
            y = y - .002
        plt.text(x, y, name, size=10,
                 horizontalalignment=horizontalalignment,
                 verticalalignment=verticalalignment,
                 bbox=dict(facecolor='w',
                           edgecolor=plt.cm.nipy_spectral(label / float(n_labels)),
                           alpha=.6))

    plt.xlim(embedding[0].min() - .15 * embedding[0].ptp(),
             embedding[0].max() + .10 * embedding[0].ptp(),)
    plt.ylim(embedding[1].min() - .03 * embedding[1].ptp(),
             embedding[1].max() + .03 * embedding[1].ptp())


    plt.savefig(os.path.join(os.path.dirname(path_to_file),'linedraw_graph.pdf'),bbox_inches="tight")

    prec_sp = gl.precision_

    save_network_graph( -prec_sp + np.diag( np.diagonal( prec_sp) ), names, os.path.join(os.path.dirname(path_to_file),"LargeNetworkNo_SP.pdf"), title="spring_sp_prec_test",layout="spring", scale= 10, node_weight=node_weights, weight = lambda x: abs(5*x)**(2.5))
    save_network_graph( gl.covariance_, names, os.path.join(os.path.dirname(path_to_file),"cov_diagram.pdf"), node_weight=node_weights, title="cov_test", scale = 8, layout = "spring" )
    save_network_graph( gl.precision_, names, os.path.join(os.path.dirname(path_to_file),"precision.pdf") , title="Precision Matrix Network", layout="spring", node_weight= node_weights)

def community_cluster(cov_sp, symbols):
	G = nx.Graph( cov_sp )
	partition = community.best_partition( G )
	for i in set(partition.values() ):
		print("Community: ",i)
		members = [ symbols[node] for node  in partition.keys() if partition[node] == i]
		print(members)


def affinity_cluster( cov_sp, symbols):
	print("Affinity cluster")
	ap = cluster.AffinityPropagation()
	ap.fit( cov_sp )
	labels = np.array(  ap.labels_ )
	for i in set(ap.labels_):
		print("Community: ",i)
		members = [ symbols[node] for node in np.nonzero( labels == i)[0]]
		print(members)
#community_cluster(gl.covariance_, top_pca_genes)
#affinity_cluster(gl.covariance_, final_gene_list)
