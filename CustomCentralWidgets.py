from CustomWidgets import *


class CentralWidgetAlignment(QWidget):
    def __init__(self, parent, raw_conc_time_graph, smoothed_conc_time_graph, ratio_dp_graph, temp_graph):
        super(self.__class__, self).__init__()
        # Add widgets
        # init the necessary contents
        self.parent = parent
        self.raw_conc_time_graph = raw_conc_time_graph
        self.smoothed_conc_time_graph = smoothed_conc_time_graph
        self.ratio_dp_graph = ratio_dp_graph
        self.temp_graph = temp_graph
        self.h_splitter_1 = QSplitter(Qt.Horizontal)
        self.h_splitter_1.addWidget(self.raw_conc_time_graph)
        self.h_splitter_1.addWidget(self.smoothed_conc_time_graph)
        self.h_splitter_2 = QSplitter(Qt.Horizontal)
        self.h_splitter_2.addWidget(self.ratio_dp_graph)
        self.h_splitter_2.addWidget(self.temp_graph)
        self.v_splitter = QSplitter(Qt.Vertical)
        self.v_splitter.addWidget(self.h_splitter_1)
        self.v_splitter.addWidget(self.h_splitter_2)
        hbox = QHBoxLayout(self)
        hbox.addWidget(self.v_splitter)
        # change some important styles
        handle = self.h_splitter_1.handle(1)
        a_layout = QVBoxLayout()
        a_layout.setSpacing(0)
        a_layout.setContentsMargins(0, 0, 0, 0)
        a_layout.addWidget(QVLine())
        handle.setLayout(a_layout)
        handle = self.h_splitter_2.handle(1)
        a_layout = QVBoxLayout()
        a_layout.setSpacing(0)
        a_layout.setContentsMargins(0, 0, 0, 0)
        a_layout.addWidget(QVLine())
        handle.setLayout(a_layout)
        handle = self.v_splitter.handle(1)
        a_layout = QHBoxLayout()
        a_layout.setSpacing(0)
        a_layout.setContentsMargins(0, 0, 0, 0)
        a_layout.addWidget(QHLine())
        handle.setLayout(a_layout)
        self.setLayout(hbox)

class CentralWidgetKappa(QWidget):
    def __init__(self, parent, kappa_graph):
        super(self.__class__, self).__init__()
        # Add widgets
        # init the necessary contents
        self.parent = parent
        self.kappa_graph = kappa_graph
        v_layout = QVBoxLayout(self)
        v_layout.addWidget(self.kappa_graph)
        self.setLayout(v_layout)
