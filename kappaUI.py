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
from timeit import default_timer as timer


class KappaInformationAndDataWidget(QWidget):
    def __init__(self, main_window=None):
        super(self.__class__, self).__init__(main_window)
        self.main_window = main_window
        self.layout = QVBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        if self.main_window.controller.is_show_all_k_points:
            self.information_table = AllKappaPointsDataTable(self.main_window)
        else:
            self.information_table = AverageKappaPointsDataTable(self.main_window)
        self.layout.addWidget(self.information_table)
        self.setLayout(self.layout)

    def resize(self, parent_width, parent_height):
        self.setFixedWidth(parent_width / 4)
        self.setFixedHeight(parent_height)
        self.information_table.resize(self.width(), self.height())

    def update_data(self):
        self.information_table.update_data()
        self.information_table.resize(self.width(), self.height())

    def show_data_for_ave_k_points(self):
        self.clear_layout(self.layout)
        self.information_table = AverageKappaPointsDataTable(self.main_window)
        self.layout.addWidget(self.information_table)

    def show_data_for_all_k_points(self):
        self.clear_layout(self.layout)
        self.information_table = AllKappaPointsDataTable(self.main_window)
        self.layout.addWidget(self.information_table)

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
        self.main_window = main_window
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def update_data(self):
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount() - 1)
        self.insertRow(self.rowCount())
        ss = AllTableHeader('ss')
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = AllTableHeader('dp')
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = AllTableHeader('K/app')
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = AllTableHeader('K/ana')
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = AllTableHeader('% devi')
        self.setCellWidget(self.rowCount() - 1, 4, devi)

        kappa_dict = self.main_window.controller.kappa_calculate_dict
        key_list = self.main_window.controller.kappa_points_data_list
        usability_list = self.main_window.controller.kappa_points_is_included_list
        for a_key in key_list:
            if kappa_dict[a_key[1]]:
                for aSS in kappa_dict[a_key[1]]:
                    if aSS[0] == a_key[0]:
                        ss = a_key[1]
                        dp = aSS[0]
                        app = aSS[1]
                        ana = aSS[2]
                        devi = aSS[3]
                        if usability_list[(dp, ss)]:
                            self.add_message(ss, dp, app, ana, devi)
                        else:
                            self.add_message(ss, dp, app, ana, devi, settings.NEGATIVE_USABILITY_COLOR)
                        break

    def add_message(self, ss, dp, app, ana, devi, color=None):
        ss = ('% .2f' % ss)
        dp = ('% .4f' % dp)
        app = ('% .4f' % app)
        ana = ('% .4f' % ana)
        devi = ('% .4f' % devi)
        self.insertRow(self.rowCount())
        ss = KappaTableItem(ss, color)
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = KappaTableItem(dp, color)
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = KappaTableItem(app, color)
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = KappaTableItem(ana, color)
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = KappaTableItem(devi, color)
        self.setCellWidget(self.rowCount() - 1, 4, devi)

    def toggle_color(self, row):
        for i in range(self.columnCount()):
            self.cellWidget(row + 1, i).toggle_color()


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
        self.main_window = main_window
        self.setAutoFillBackground(True)
        palette = QPalette()
        self.SelectItems
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        palette.setColor(QPalette.Highlight, settings.NEGATIVE_USABILITY_COLOR)
        palette.setColor(QPalette.HighlightedText, settings.NEGATIVE_USABILITY_COLOR)
        self.setPalette(palette)

    def update_data(self):
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount() - 1)
        self.insertRow(self.rowCount())
        ss = AllTableHeader('ss')
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = AllTableHeader('dp')
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = AllTableHeader('K/app')
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = AllTableHeader('K/ana')
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = AllTableHeader('% devi')
        self.setCellWidget(self.rowCount() - 1, 4, devi)

        data_dict = self.main_window.controller.alpha_pinene_dict
        key_list = self.main_window.controller.kappa_points_data_list
        for a_key in key_list:
            a_key = a_key[1]
            if data_dict[a_key]:
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
        ss = KappaTableItem(ss)
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = KappaTableItem(dp)
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = KappaTableItem(app)
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = KappaTableItem(ana)
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = KappaTableItem(devi)
        self.setCellWidget(self.rowCount() - 1, 4, devi)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height)
        self.setFixedWidth(parent_width)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def reset_data(self):
        for i in range(self.rowCount()):
            for j in range(self.columnCount()):
                self.cellWidget(i, j).toggle_color()

    def toggle_color(self, row):
        for i in range(self.columnCount()):
            self.cellWidget(row + 1, i).toggle_color()


class KappaGraphWidget(QWidget):
    def resize(self, parent_width, parent_height):
        self.setFixedWidth(parent_width * 3 / 4)
        self.setFixedHeight(parent_height)
        self.kappa_graph_view.resize(self.width(), self.height())
        self.control_widget.resize(self.width(), self.height())

    def __init__(self, main_window=None):
        super(self.__class__, self).__init__(main_window)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.kappa_graph_view = KappaFigureCanvas()
        self.control_widget = KappaControlTabWidget(main_window)
        self.layout.addWidget(self.kappa_graph_view)
        self.layout.addWidget(self.control_widget)

    def update_figure(self, figure):
        self.kappa_graph_view.update_figure(figure)


class KappaControlTabWidget(QWidget):
    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 1 / 10)
        self.setFixedWidth(parent_width)
        self.show_parameters_button.resize(self.width(), self.height())
        self.toggle_average_all_k_points_button.resize(self.width(), self.height())
        self.toggle_k_point_status_button.resize(self.width(),self.height())
        self.export_data_button.resize(self.width(),self.height())

    def __init__(self, main_window=None):
        self.main_window = main_window
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0, 10, 60, 10)

        self.show_parameters_button = CustomButton("Parameters", main_window)
        self.toggle_k_point_status_button = CustomButton("Disable", main_window)
        self.toggle_average_all_k_points_button = CustomButton("Ave Points", main_window)
        self.export_data_button = CustomButton("Export To CSV", main_window)

        self.show_parameters_button.clicked.connect(self.on_click_show_parameters)
        self.toggle_average_all_k_points_button.clicked.connect(self.on_click_toggle_all_average_k_points)
        self.toggle_k_point_status_button.clicked.connect(self.on_click_toggle_k_point_status_button)
        self.export_data_button.clicked.connect(self.on_click_export_data_button)

        self.layout.addWidget(self.show_parameters_button)
        self.layout.addWidget(self.toggle_k_point_status_button)
        self.layout.addWidget(self.toggle_average_all_k_points_button)
        self.layout.addWidget(self.export_data_button)

        self.setLayout(self.layout)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, settings.controlAreaBackgroundColor)
        self.setPalette(palette)

    def on_click_show_parameters(self):
        controller = self.main_window.controller
        kappaVars = (
            controller.sigma, controller.temp, controller.dd, controller.iKappa, controller.dd2, controller.iKappa2,
            controller.solubility)
        sig = kappaVars[0]
        temp = kappaVars[1]
        dd1 = kappaVars[2]
        iKappa1 = kappaVars[3]
        dd2 = kappaVars[4]
        iKappa2 = kappaVars[5]
        solu = kappaVars[6]
        kappa_var_dialog = KappaVarDialog(sig, temp, dd1, iKappa1, dd2, iKappa2, solu, self.main_window)
        kappa_var_dialog.exec_()

    def on_click_toggle_all_average_k_points(self):
        if (self.main_window.controller.is_show_all_k_points):
            self.main_window.controller.is_show_all_k_points = False
            self.toggle_average_all_k_points_button.setText("All Points")
            self.main_window.centralWidget().info_widget.show_data_for_ave_k_points()
            self.toggle_k_point_status_button.hide()
        else:
            self.main_window.controller.is_show_all_k_points = True
            self.toggle_average_all_k_points_button.setText("Ave Points")
            self.main_window.centralWidget().info_widget.show_data_for_all_k_points()
            self.toggle_k_point_status_button.show()
        # reset certain attributes of kappa to draw correctly
        self.main_window.controller.current_kappa_point_index = None
        self.main_window.controller.update_kappa_info_and_graph()

    def on_click_toggle_k_point_status_button(self):
        self.main_window.controller.toggle_exclude_include_kappa_point()

    def on_click_export_data_button(self):
        self.main_window.controller.export_to_csv()


class KappaFigureCanvas(FigureCanvas):
    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 9 / 10 + 5)
        self.setFixedWidth(parent_width)

    def __init__(self, main_window=None):
        fig = Figure(facecolor=settings.GRAPH_BACKGROUND_COLOR)
        super(self.__class__, self).__init__(fig)

    def update_figure(self, new_figure):
        if new_figure is None:
            return
        if self.figure != new_figure:
            self.figure = new_figure
            h = self.height()
            self.setFixedHeight(h / 2)
            self.setFixedHeight(h)
        else:
            self.draw()
            self.flush_events()


class KappaVarDialog(QDialog):
    def __init__(self, sigma, temp, dd1, iKappa1, dd2, iKapp2, solu, main_window=None):
        super(self.__class__, self).__init__()
        self.mainWindow = main_window
        self.formLayout = QFormLayout()
        self.sigmaLine = QLabel(str(sigma))
        self.tempLine = QLabel(str(temp))
        self.dd1Line = QLabel(str(dd1))
        self.iKappa1Line = QLabel(str(iKappa1))
        self.dd2Line = QLabel(str(dd2))
        self.iKappa2Line = QLabel(str(iKapp2))
        self.soluLine = QLabel(str(solu))
        self.formLayout.addRow(self.tr("&Sigma"), self.sigmaLine)
        self.formLayout.addRow(self.tr("&Temperature"), self.tempLine)
        self.formLayout.addRow(self.tr("&dry diameter(1)"), self.dd1Line)
        self.formLayout.addRow(self.tr("&iKappa(1)"), self.iKappa1Line)
        self.formLayout.addRow(self.tr("&dry diameter(2)"), self.dd2Line)
        self.formLayout.addRow(self.tr("&iKappa(2)"), self.iKappa2Line)
        self.formLayout.addRow(self.tr("&solubility"), self.soluLine)
        self.setLayout(self.formLayout)
