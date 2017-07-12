from InfoTableStyle import *
from CustomButton import *
import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time


class KappaInformationAndDataWidget(QWidget):
    def __init__(self, main_window=None):
        super(self.__class__, self).__init__(main_window)
        self.main_window = main_window
        self.layout = QVBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.current_data_table = KappaParamsDataTable(main_window)
        self.layout.addWidget(self.current_data_table)
        self.setLayout(self.layout)

    def resize(self, parent_width, parent_height):
        self.setFixedWidth(parent_width / 4)
        self.setFixedHeight(parent_height)
        self.current_data_table.resize(self.width(), self.height())

    def update_data(self):
        self.current_data_table.update_data()
        self.current_data_table.resize(self.width(), self.height())

    def show_data_for_ave_k_points(self):
        self.clear_layout(self.layout)
        self.current_data_table = AverageKappaPointsDataTable(self.main_window)
        self.layout.addWidget(self.current_data_table)
        self.update_data()

    def show_data_for_all_k_points(self):
        self.clear_layout(self.layout)
        self.current_data_table = AllKappaPointsDataTable(self.main_window)
        self.layout.addWidget(self.current_data_table)
        self.update_data()

    def change_to_parameters_data_table(self):
        self.clear_layout(self.layout)
        self.current_data_table = KappaParamsDataTable(self.main_window)
        self.layout.addWidget(self.current_data_table)
        self.update_data()

    def clear_layout(self, layout):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clear_layout(item.layout())


class AllKappaPointsDataTable(QTableWidget):
    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height)
        self.setFixedWidth(parent_width)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def __init__(self, main_window=None):
        QTableWidget.__init__(self)
        self.setColumnCount(5)
        self.setRowCount(0)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.setFrameStyle(QFrame.NoFrame)
        self.mainWindow = main_window
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def update_data(self):
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount()-1)
        self.insertRow(self.rowCount())
        ss = TableHeader('ss')
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = TableHeader('dp')
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = TableHeader('K/app')
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = TableHeader('K/ana')
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = TableHeader('% devi')
        self.setCellWidget(self.rowCount() - 1, 4, devi)

        kappa_dict = self.mainWindow.controller.kappa_calculate_dict
        for a_key in kappa_dict.keys():
            if kappa_dict[a_key]:
                for aSS in kappa_dict[a_key]:
                    ss = a_key
                    dp = aSS[0]
                    app = aSS[1]
                    ana = aSS[2]
                    devi = aSS[3]
                    self.add_message(ss, dp, app, ana, devi)

    def add_message(self, ss, dp, app, ana, devi):
        ss = ('% .2f' % ss)
        dp = ('% .4f' % dp)
        app = ('% .4f' % app)
        ana = ('% .4f' % ana)
        devi = ('% .4f' % devi)
        self.insertRow(self.rowCount())
        ss = SingleTableItem(ss)
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = SingleTableItem(dp)
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = SingleTableItem(app)
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = SingleTableItem(ana)
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = SingleTableItem(devi)
        self.setCellWidget(self.rowCount() - 1, 4, devi)


class AverageKappaPointsDataTable(QTableWidget):
    def __init__(self, main_window=None):
        QTableWidget.__init__(self)
        self.setColumnCount(5)
        self.setRowCount(0)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.setFrameStyle(QFrame.NoFrame)
        self.mainWindow = main_window
        self.setAutoFillBackground(True)
        palette = QPalette()
        self.SelectItems
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        palette.setColor(QPalette.Highlight, settings.NEGATIVE_USABILITY_BUTTON_COLOR)
        palette.setColor(QPalette.HighlightedText,settings.NEGATIVE_USABILITY_BUTTON_COLOR)
        self.setPalette(palette)

    def update_data(self):
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount()-1)
        self.insertRow(self.rowCount())
        ss = TableHeader('ss')
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = TableHeader('dp')
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = TableHeader('K/app')
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = TableHeader('K/ana')
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = TableHeader('% devi')
        self.setCellWidget(self.rowCount() - 1, 4, devi)

        data_dict = self.mainWindow.controller.alpha_pinene_dict
        for a_key in data_dict.keys():
            ss = a_key
            dp = data_dict[a_key][0]
            app = data_dict[a_key][2]
            ana = data_dict[a_key][4]
            devi = data_dict[a_key][6]
            self.add_message(ss, dp, app, ana, devi)

    def add_message(self, ss, dp, app, ana, devi):
        ss = ('% .2f' % ss)
        dp = ('% .4f' % dp)
        app = ('% .4f' % app)
        ana = ('% .4f' % ana)
        devi = ('% .4f' % devi)
        self.insertRow(self.rowCount())
        ss = SingleTableItem(ss)
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = SingleTableItem(dp)
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = SingleTableItem(app)
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = SingleTableItem(ana)
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = SingleTableItem(devi)
        self.setCellWidget(self.rowCount() - 1, 4, devi)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height)
        self.setFixedWidth(parent_width)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.cellWidget(2,2).toggle_color()




class KappaParamsDataTable(QTableWidget):
    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height)
        self.setFixedWidth(parent_width)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def __init__(self, main_window=None):
        QTableWidget.__init__(self)
        self.setColumnCount(1)
        self.setRowCount(0)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.resizeRowsToContents()
        self.setFrameStyle(QFrame.NoFrame)
        self.mainWindow = main_window

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def update_data(self):
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount() - 1)

        header = TableHeader("Variables")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        controller = self.mainWindow.controller
        kappaVars = (controller.sigma, controller.temp, controller.dd, controller.iKappa, controller.dd2, controller.iKappa2, controller.solubility)
        s = kappaVars[0]
        t = kappaVars[1]
        d1 = kappaVars[2]
        i1 = kappaVars[3]
        d2 = kappaVars[4]
        i2 = kappaVars[5]
        solu = kappaVars[6]
        # TODO: fix this
        self.add_message("Sigma", s)
        self.add_message("Temperature (k)", t)
        self.add_message("dry diamater(1) (nm)", d1)
        self.add_message("iKppa(1)", i1)
        self.add_message("dry diamater(2) (nm)", d2)
        self.add_message("iKppa(2)", i2)
        self.add_message("Solubility", solu)

    def add_message(self, field, message):
        message = str(message)
        item = TableItem(field, message)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount()-1, 0, item)


class KappaGraphWidget(QWidget):
    def resize(self, parent_width, parent_height):
        self.setFixedWidth(parent_width * 3 / 4)
        self.setFixedHeight(parent_height)
        self.graphView.resize(self.width(), self.height())
        self.controlWidget.resize(self.width(), self.height())

    def __init__(self, main_window=None):
        super(self.__class__, self).__init__(main_window)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.graphView = KappaFigureCanvas()
        self.controlWidget = KappaControlTabWidget(main_window)
        self.layout.addWidget(self.graphView)
        self.layout.addWidget(self.controlWidget)


class KappaControlTabWidget(QWidget):
    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 1 / 10)
        self.setFixedWidth(parent_width)
        self.show_parameters_button.resize(self.width(), self.height())
        self.toggle_average_all_k_points_button.resize(self.width(), self.height())

    def __init__(self, main_window=None):
        self.mainWindow = main_window
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0, 10, 60, 10)

        self.show_parameters_button = CustomButton("Parameters", main_window)
        self.toggle_average_all_k_points_button = CustomButton("Ave Points", main_window)

        self.show_parameters_button.clicked.connect(self.on_click_show_parameters)
        self.toggle_average_all_k_points_button.clicked.connect(self.on_click_toggle_all_average_k_points)

        self.layout.addWidget(self.show_parameters_button)
        self.layout.addWidget(self.toggle_average_all_k_points_button)

        self.setLayout(self.layout)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, settings.controlAreaBackgroundColor)
        self.setPalette(palette)
        # self.toggle_all_less_k_lines_button = CustomButton("Less Lines", main_window)
        # self.show_raw_data_button = CustomButton("Kappa Data", main_window)
        # self.show_graph_data_button = CustomButton("Graph Data", main_window)
        # self.show_raw_data_button.clicked.connect(self.on_click_show_raw_data)
        # self.show_graph_data_button.clicked.connect(self.on_click_show_graph)
        # self.layout.addWidget(self.toggle_all_less_k_lines_button)
        # self.layout.addWidget(self.show_raw_data_button)
        # self.layout.addWidget(self.show_graph_data_button)
        # self.toggle_all_less_k_lines_button.clicked.connect(self.on_click_toggle_k_lines)

    def on_click_show_parameters(self):
        self.mainWindow.centralWidget().info_widget.change_to_parameters_data_table()

    def on_click_toggle_all_average_k_points(self):
        if (self.mainWindow.controller.is_show_all_k_points):
            self.mainWindow.controller.is_show_all_k_points = False
            self.toggle_average_all_k_points_button.setText("All Points")
            self.mainWindow.centralWidget().info_widget.show_data_for_ave_k_points()
        else:
            self.mainWindow.centralWidget().info_widget.show_data_for_all_k_points()
            self.mainWindow.controller.is_show_all_k_points = True
            self.toggle_average_all_k_points_button.setText("Ave Points")
        self.mainWindow.controller.draw_kappa_graph()

    # def on_click_toggle_k_lines(self):
    #     if (self.main_window.controller.is_show_all_k_lines):
    #         self.main_window.controller.is_show_all_k_lines = False
    #         self.toggle_all_less_k_lines_button.setText("All Lines")
    #     else:
    #         self.main_window.controller.is_show_all_k_lines = True
    #         self.toggle_all_less_k_lines_button.setText("Less Lines")
    #     self.main_window.controller.draw_kappa_graph()

    # def on_click_show_graph(self):
    #     self.main_window.centralWidget().info_widget.show_data_for_ave_k_points()

    # def on_click_show_raw_data(self):
    #     self.main_window.centralWidget().info_widget.show_data_for_all_k_points()


class KappaFigureCanvas(FigureCanvas):
    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 9 / 10 + 5)
        self.setFixedWidth(parent_width)

    def __init__(self, main_window=None):
        fig = Figure(facecolor=settings.GRAPH_BACKGROUND_COLOR)
        super(self.__class__, self).__init__(fig)

    def update_figure(self, figure):
        self.figure = figure
        self.draw()

        # A hack to make the figure update to the size of the Figure Canvas
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)
