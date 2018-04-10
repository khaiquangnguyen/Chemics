from MainView import *
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize as opt
import scipy.stats
import scipy.signal
from Exceptions import *
import gc
from HelperFunctions import *
import FastDpCalculator
import settings as CONST
import matplotlib.patches as mpatches
import os
from PySide.QtCore import SIGNAL
from PySide.QtCore import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
matplotlib.style.use('ggplot')

class ConcOverTimeRawDataGraph(FigureCanvas):
    def __init__(self, parent = None):
        self.fig, self.ax = plt.subplots(facecolor=CONST.GRAPH_BACKGROUND_COLOR)
        super(self.__class__, self).__init__(self.fig)
        self.ax.set_autoscale_on(True)
        self.ax.axes.set_frame_on(False)
        self.ax.grid(color=CONST.GRID_COLOR)
        self.ax.axhline(0, color=CONST.AX_LINE_COLOR, linewidth=4)
        self.ax.axvline(0, color=CONST.AX_LINE_COLOR, linewidth=4)
        self.ax.tick_params(color=CONST.AX_TICK_COLOR, which='both', labelcolor=CONST.AX_TICK_COLOR, labelsize=CONST.AX_TICK_SIZE)
        self.ax.set_xlabel("Scan time(s)", color=CONST.LABEL_COLOR, size=CONST.LABEL_SIZE)
        self.ax.set_ylabel("Concentration (1/cm3)", color=CONST.LABEL_COLOR, size=CONST.LABEL_SIZE)
        self.smps_points, = self.ax.plot([],[], linewidth=4, color=CONST.SMPS_LINE_COLOR, label="Raw SMPS")
        self.ccnc_points, = self.ax.plot([],[],linewidth=4, color = CONST.CCNC_LINE_COLOR, label = "Raw CCNC")
        handles, labels = self.ax.get_legend_handles_labels()
        legend = self.ax.legend(handles, labels, loc=2, fontsize=CONST.LEGEND_FONT_SIZE)
        legend.get_frame().set_alpha(0.3)
        self.ax.set_title("Raw SMPS and CCNC concentration over time", color=CONST.TITLE_COLOR, size=CONST.TITLE_SIZE)

    def update_graph(self, a_scan, smooth_ccnc = False):
        smps_counts = a_scan.raw_smps_counts
        ccnc_counts = a_scan.raw_ccnc_counts
        # if this graph is to show the ccnc data with smoothing
        # mostly use for showcase, who knows what the usage are
        if smooth_ccnc:
            ccnc_counts = smooth(ccnc_counts)
        # get the number of data points
        num_data_pts = max(len(smps_counts),len(ccnc_counts))
        # make up for the lost data points
        smps_counts = fill_zeros_to_end(smps_counts,num_data_pts)
        ccnc_counts = fill_zeros_to_end(ccnc_counts,num_data_pts)
        smps_counts = smps_counts[:num_data_pts]
        ccnc_counts = ccnc_counts[:num_data_pts]
        # find the limits of y-axis
        min_y = numpy.amin([smps_counts, ccnc_counts])
        max_y = numpy.amax([smps_counts, ccnc_counts])
        # set the limits
        self.ax.axes.set_ylim([min_y, max_y])
        self.ax.axes.set_xlim([0, num_data_pts])
        # set new data
        x_axis = numpy.arange(num_data_pts)
        self.smps_points.set_xdata(x_axis)
        self.smps_points.set_ydata(smps_counts)
        self.ccnc_points.set_xdata(x_axis)
        self.ccnc_points.set_ydata(ccnc_counts)
        self.draw()
        self.flush_events()


class ConcOverTimeSmoothGraph(FigureCanvas):
    def __init__(self, parent = None):
        self.fig, self.ax = plt.subplots(facecolor=CONST.GRAPH_BACKGROUND_COLOR)
        super(self.__class__, self).__init__(self.fig)
        self.ax.set_autoscale_on(True)
        self.ax.axes.set_frame_on(False)
        self.ax.grid(color=CONST.GRID_COLOR)
        self.ax.axhline(0, color=CONST.AX_LINE_COLOR, linewidth=4)
        self.ax.axvline(0, color=CONST.AX_LINE_COLOR, linewidth=4)
        self.ax.tick_params(color=CONST.AX_TICK_COLOR, which='both', labelcolor=CONST.AX_TICK_COLOR, labelsize=CONST.AX_TICK_SIZE)
        self.ax.set_xlabel("Scan time(s)", color=CONST.LABEL_COLOR, size=CONST.LABEL_SIZE)
        self.ax.set_ylabel("Concentration (1/cm3)", color=CONST.LABEL_COLOR, size=CONST.LABEL_SIZE)
        self.smps_points, = self.ax.plot([],[], linewidth=4, color=CONST.SMPS_LINE_COLOR, label="SMPS")
        self.ccnc_points, = self.ax.plot([],[],linewidth=4, color = CONST.CCNC_LINE_COLOR, label = "CCNC")
        handles, labels = self.ax.get_legend_handles_labels()
        legend = self.ax.legend(handles, labels, loc=2, fontsize=CONST.LEGEND_FONT_SIZE)
        legend.get_frame().set_alpha(0.3)
        self.ax.set_title("SMPS and CCNC concentration over scan time", color=CONST.TITLE_COLOR, size=CONST.TITLE_SIZE)

    def update_graph(self, a_scan):
        smps_counts = a_scan.processed_smps_counts
        ccnc_counts = a_scan.processed_ccnc_counts
        # get the number of data points
        num_data_pts = a_scan.duration
        # make up for the lost data points
        smps_counts = fill_zeros_to_end(smps_counts, num_data_pts)
        ccnc_counts = fill_zeros_to_end(ccnc_counts, num_data_pts)
        smps_counts = smps_counts[:num_data_pts]
        ccnc_counts = ccnc_counts[:num_data_pts]
        ccnc_counts = smooth(ccnc_counts)
        # find the limits of y-axis
        min_y = numpy.amin([smps_counts, ccnc_counts])
        max_y = numpy.amax([smps_counts, ccnc_counts])
        # set the limits
        self.ax.axes.set_ylim([min_y, max_y])
        self.ax.axes.set_xlim([0, num_data_pts])
        # set new data
        x_axis = numpy.arange(num_data_pts)
        self.smps_points.set_xdata(x_axis)
        self.smps_points.set_ydata(smps_counts)
        self.ccnc_points.set_xdata(x_axis)
        self.ccnc_points.set_ydata(ccnc_counts)
        self.draw()
        self.flush_events()


class RatioOverDiameterGraph(FigureCanvas):
    def __init__(self,parent=None):
        self.fig, self.ax = plt.subplots(facecolor=CONST.GRAPH_BACKGROUND_COLOR)
        super(self.__class__, self).__init__(self.fig)
        self.ax.axes.set_frame_on(False)
        self.ax.grid(color='0.5')
        self.ax.axhline(0, color=CONST.AX_LINE_COLOR, linewidth=2)
        self.ax.axvline(0, color=CONST.AX_LINE_COLOR, linewidth=4)
        self.ax.axhline(1, color=str(float(CONST.AX_LINE_COLOR) + 0.1), linewidth=2, linestyle='dashed')
        self.ax.tick_params(color=CONST.AX_TICK_COLOR, which='both', labelcolor=CONST.AX_TICK_COLOR, labelsize=CONST.AX_TICK_SIZE)
        self.ax.set_xlabel("Diameter (nm)", color=CONST.LABEL_COLOR, size=CONST.LABEL_SIZE)
        self.ax.set_ylabel("CCNC/SMPS", color=CONST.LABEL_COLOR, size=CONST.LABEL_SIZE)
        self.ax.set_title("CCNC/SMPS over Dry Diameter", color=CONST.TITLE_COLOR, size=CONST.TITLE_SIZE)
        self.ccn_cn_ratio_points, = self.ax.plot([0],[0], 'o',
                                            color=CONST.CCNC_SMPS_POINT_COLOR, mew=0.5,
                                            ms=9, label="CCNC/SMPS")
        handles, labels = self.ax.get_legend_handles_labels()
        legend = self.ax.legend(handles, labels, loc=5, fontsize=CONST.LEGEND_FONT_SIZE)
        legend.get_frame().set_alpha(0.3)

    def update_graph(self,a_scan):
        # extract the data
        ccnc = a_scan.processed_ccnc_counts
        smps = a_scan.processed_smps_counts
        ave_smps_dp = a_scan.ave_smps_diameters
        # calculate the ratio
        ratio = []
        for i in range(len(smps)):
            ratio.append(safe_div(ccnc[i],smps[i]))
        # we only need data within 0 to 1.5. Anything beyond that should be errors, so we don't need to graph them
        # set the limits of the axes
        y_max = numpy.amin([1.5, numpy.amax(ratio)])
        self.ax.axes.set_ylim([0, y_max])
        x_min = numpy.amin(ave_smps_dp)
        x_max = numpy.amax(ave_smps_dp)
        self.ax.axes.set_xlim([x_min, x_max])
        self.ccn_cn_ratio_points.set_xdata(ave_smps_dp)
        self.ccn_cn_ratio_points.set_ydata(ratio)
        self.draw()
        self.flush_events()


class TemperatureGraph(FigureCanvas):
    def __init__(self,parent=None):
        self.fig, self.ax = plt.subplots(facecolor=CONST.GRAPH_BACKGROUND_COLOR)
        # set up the figure and axes
        self.fig, self.ax = plt.subplots(facecolor=CONST.GRAPH_BACKGROUND_COLOR)
        super(self.__class__, self).__init__(self.fig)
        self.ax.set_autoscale_on(True)
        self.ax.set_frame_on(False)
        self.ax.grid(color=CONST.GRID_COLOR)
        self.ax.axhline(0, color=CONST.AX_LINE_COLOR, linewidth=4)
        self.ax.axvline(0, color=CONST.AX_LINE_COLOR, linewidth=4)
        self.ax.tick_params(color=CONST.AX_TICK_COLOR, which='both', labelcolor=CONST.AX_TICK_COLOR, labelsize=CONST.AX_TICK_SIZE)
        self.ax.set_xlabel("Scan time(s)", color=CONST.LABEL_COLOR, size=CONST.LABEL_SIZE)
        self.ax.set_ylabel("Temperature", color=CONST.LABEL_COLOR, size=CONST.LABEL_SIZE)
        self.ax.set_title("Temperature over scan time", color=CONST.TITLE_COLOR, size=CONST.TITLE_SIZE)
        # set up empty data lines
        self.t1_line, = self.ax.plot([0], [0], linewidth=5, color=CONST.T1_LINE_COLOR, label="T1")
        self.t2_line, = self.ax.plot([0], [0], linewidth=5, color=CONST.T2_LINE_COLOR, label="T2")
        self.t3_line, = self.ax.plot([0], [0], linewidth=5, color=CONST.T3_LINE_COLOR, label="T3")
        handles, labels = self.ax.get_legend_handles_labels()
        legend = self.ax.legend(handles, labels, fontsize=CONST.LEGEND_FONT_SIZE)
        legend.get_frame().set_alpha(0.3)

    def update_graph(self, a_scan):
        """
        Update the temperature graph based on the data of a scan.
        :param a_scan: OneScan object. Contains temperatures information
        """
        t1s = a_scan.processed_T1s
        t2s = a_scan.processed_T2s
        t3s = a_scan.processed_T3s
        # get the number of data points
        num_data_pts = a_scan.duration
        t1s = fill_zeros_to_end(t1s,num_data_pts)
        t2s = fill_zeros_to_end(t2s,num_data_pts)
        t3s = fill_zeros_to_end(t3s,num_data_pts)
        # find the limits of y-axis
        min_y = numpy.amin([t1s,t2s,t3s])
        max_y = numpy.amax([t1s,t2s,t3s])
        # widen the range of data a bit
        min_y = min_y - max_y/10
        max_y = max_y + max_y / 10
        # set the limits
        self.ax.axes.set_ylim([min_y,max_y])
        self.ax.axes.set_xlim([0, num_data_pts])
        # set new data
        x_axis = numpy.arange(num_data_pts)
        self.t1_line.set_xdata(x_axis)
        self.t1_line.set_ydata(t1s)
        self.t2_line.set_xdata(x_axis)
        self.t2_line.set_ydata(t2s)
        self.t3_line.set_xdata(x_axis)
        self.t3_line.set_ydata(t3s)
        self.draw()
        self.flush_events()




#
# class SigmoidFitGraph:
#
# class ScanSummaryGraph:
#
# class KappaGraph:
