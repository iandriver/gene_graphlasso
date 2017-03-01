#!/usr/bin/python3
from lasso_functions import *
from qt5_gui_gl import Graph_UI
from PyQt5.QtWidgets import QApplication




def main():
    app = QApplication(sys.argv)
    qt5_data = Graph_UI()
    app.exec_()
    top_matrix, top_matrix_cell, final_gene_list = make_pcor_network(qt5_data)
    gl, embedding, names = make_graphlasso(top_matrix, final_gene_list)
    node_weights = find_node_weight(top_matrix, final_gene_list)
    plot_networks(gl, embedding, names, node_weights, qt5_data)


if __name__ == '__main__':
    main()
