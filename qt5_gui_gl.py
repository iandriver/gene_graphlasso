import sys
from PyQt5.QtWidgets import (QWidget,QMainWindow, QTextEdit,
    QAction, QFileDialog, QApplication, QPushButton, QLineEdit, QGridLayout, QCheckBox, QInputDialog, QLabel)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt, QCoreApplication


class Graph_UI(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):

        grid = QGridLayout()
        grid.setSpacing(10)

        self.kmeans = False
        self.btn = QPushButton('Select File', self)
        #self.btn.move(30,50)
        self.btn.clicked.connect(self.fileDialog)



        self.filepath = QLineEdit(self)
        #self.le.move(20, 10)
        self.filepath.setFixedWidth(300)
        self.filepath.setText('')

        self.iter_label = QLabel("# of search iterations: ")
        self.iterations = QLineEdit(self)
        #self.le.move(20, 10)
        self.iterations.setFixedWidth(50)
        self.iterations.setText('1')

        self.qbtn = QPushButton('Run', self)
        self.qbtn.clicked.connect(self.get_itertext)
        self.qbtn.clicked.connect(QCoreApplication.instance().quit)
        self.qbtn.resize(self.qbtn.sizeHint())
        self.qbtn.move(50, 50)

        grid.addWidget(self.filepath, 0,0)
        grid.addWidget(self.btn, 1,0)
        grid.addWidget(self.iter_label, 2,0)
        grid.addWidget(self.iterations, 2,1)
        grid.addWidget(self.qbtn, 3,0)

        self.setLayout(grid)

        #self.setGeometry(300, 300, 350, 300)
        self.setWindowTitle('Gene Graph')
        self.show()


    def fileDialog(self):

        fname = QFileDialog.getOpenFileName(self, 'Open file', '.')

        if fname[0]:
            self.filepath.setText(fname[0])
        print(self.filepath.text())
    def get_itertext(self):
        text = self.iterations.text()
        self.iterations.setText(text)
        print(text)
