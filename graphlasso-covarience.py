#!/usr/bin/python3
from lasso_functions import *
from qt5_gui_gl import Graph_UI
from PyQt5.QtWidgets import QApplication




def main():
    app = QApplication(sys.argv)
    qt5_data = Graph_UI()
    app.exec_()
    log2_df_cell, num_iterations, gene_index_list = read_qt5_input(qt5_data)
    top_matrix_gene, top_matrix_cell, final_gene_list = find_top_genes_from_pcor_network(log2_df_cell, num_iterations, gene_index_list)
    gl_gene, embedding_gene, names_gene = make_graphlasso(top_matrix_gene, final_gene_list)
    gl_cell, embedding_cell, names_cell = make_graphlasso(top_matrix_cell, top_matrix_cell.columns.tolist())
    node_weights_gene = find_node_weight(top_matrix_gene, final_gene_list)
    node_weights_cell = find_node_weight(top_matrix_cell, top_matrix_cell.columns.tolist())
    plot_networks(gl_gene, embedding_gene, names_gene, node_weights_gene, qt5_data)
    plot_networks(gl_cell, embedding_cell, names_cell, node_weights_cell, qt5_data, gene_matrix=False)
    community_cluster(top_matrix_cell,qt5_data)
    affinity_cluster(top_matrix_cell,qt5_data)


if __name__ == '__main__':
    main()
