"""
------------------------------------------------------------------------------------------------------------------------
Creates figures for manuscript
------------------------------------------------------------------------------------------------------------------------
Creates geometry and morphological characteristic plot
Creates plots for comparing comsol vs simgraph
- pressure plot -> xlabel: segment position from the origin ylabel: pressure value
- volummetric flow plot -> xlabel: segment position from the origin ylabel: velocity value
"""
import os
import string
import vtk
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib
import math
import cv2

from typing import Dict, Any, List
from pathlib import Path
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from image_processing.centerline.read_mesh import mesh_network
from matpancreas.utils import io_utils_py, graph_utils_py
from types import SimpleNamespace
from vedo import * #Lines, Points, show, Plotter, interactive, settings, Picture, screenshot, buildLUT, Axes, Text
from vedo.pyplot import plot
from matpancreas.pancreas.analysis.compare_lines import get_static_data, get_dynamic_data, get_bounds_static
from matpancreas.pancreas.analysis.compare_lines_pscan import get_static_data as get_pscan_static_data,\
    get_dynamic_data as get_pscan_dynamic_data
from matpancreas.settings_model import RESULTS_SIMGRAPH_FILE, \
    RESULTS_COMSOL_FILE, \
    comsol, \
    RESULTS_PLOTS_DIR, \
    RESULTS_PSCAN_FILE, \
    test_case as test, \
    ESPECIES, \
    single, \
    pscan, \
    pressure_scan, \
    RESULTS_SIMGRAPH_DIR, \
    RESULTS_COMSOL_DIR, \
    test_cases, pb, nspecies, SPECIES, RESULTS_GSCAN_FILE

from matpancreas.utils.graph_utils_py import get_graph, \
    create_graph, \
    get_eucledian_lengths_wrt_origin, \
    draw_graph3d_vedo, \
    plot_histogram, convert_graph_to_df, get_vals_segment

# vedo settings
from matpancreas.utils.io_utils_py import write_output

settings.screenshotTransparentBackground = True
settings.screeshotScale = 1
settings.screeshotLargeImage = True
settings.showRendererFrame = False

# global settings for plots
# plt.rcParams.update({
#         'axes.labelsize': 'large',
#         'axes.labelweight': 'bold'
#     })

# -------------- test label settings ---------------------------------------------------------------------------
left, width = .35, .5
bottom, height = .40, .5
right = left + width
top = bottom + height
fig_labels = dict(zip(test_cases, ["c", "a", "b", "d"]))
# --------------------------------------------------------------------------------------------------------------


def computeTicks (x, step=5):
    """
    Computes domain with given step encompassing series x
    @ params
    x    - Required - A list-like object of integers or floats
    step - Optional - Tick frequency
    Reference:
    https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib
    """
    xMax, xMin = math.ceil(max(x)), math.floor(min(x))
    dMax, dMin = xMax + abs((xMax % step) - step) + (step if (xMax % step != 0) else 0), xMin - abs((xMin % step))
    return range(dMin, dMax, step)


def plot_histogram(ax, data: List, binwidth: int, density=False, bin=False):
    if bin:
        ax.hist(
            data,
            density=density,
            color='k',
            bins=np.arange(min(data), max(data) + binwidth, binwidth)
        )
    else:
        ax.hist(
            data,
            density=density,
            color='k'
        )


def draw_weighted_connectivity_matrix(weights: List[List], test_cases):
    """
    This function creates a heatmap
    - incidence matrix, M (V x E), is created
    - w is weights of edges
    - W  = diag(w)
    - M x W -> column scaled incidence matirx
    : param G , graph
    :return:
    ref: https://github.com/mne-tools/mne-python/issues/5693 (to assign white color to zero)
    https://stackoverflow.com/questions/55306803/matplotlib-colorbar-need-to-force-scientific-notation-with-exponent-at-top-of-b
    """

    fig, axes = plt.subplots(nrows=1, ncols=len(test_cases), figsize=(16, 5))
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    for test_case, weight in zip(test_cases, weights):
        i = test_cases.index(test_case)
        sheet = test_case.split("_pb")[0]

        output = get_graph(sheet=sheet)
        ed_ls = output['ed_ls']
        G = create_graph(output=output)
        print(G.nodes())
        M = nx.incidence_matrix(G, nodelist=sorted(G.nodes), edgelist=ed_ls, oriented=True)
        W = np.diag(weight)
        MW = M @ W
        vmin, vmax = np.min(MW), np.max(MW)

        cbformat = matplotlib.ticker.ScalarFormatter()  # create the formatter
        cbformat.set_powerlimits((-2, 2))  # set the limits for sci. not.
        im = axes[i].matshow(
                    MW, cmap=plt.cm.bwr, #bwr_r
                    vmin=-abs(max(vmin, vmax)),
                    vmax=abs(max(vmin, vmax))
                    )
        axes[i].set_ylabel('Nodes')
        axes[i].set_xlabel('Edges')
        axes[i].text(
            0.9, 0.9,
            '(' + fig_labels[test_case] + ')',
            # string.ascii_uppercase[idx],
            transform=axes[i].transAxes,
            # weight='bold',
            size=14
        )
        plt.colorbar(
                     im,
                     ax=axes[i],
                     fraction=0.057,
                     pad=0.06,
                     format=cbformat, #'%.0e'
                     orientation="horizontal"
                     ).set_label('qdot ($\mu$m$^3$/s)', rotation=0)

    return fig


def get_plot_data(fpath_s: Path, fpath_c: Path, measure: str, sort_by_distance=False) -> Dict:
    """
    :param sort_by_distance:
    :param fpath_s:
    :param fpath_c:
    :param measure:
    :return:
    """
    labels = {'pressure': {'xlabel': 'Node', 'ylabel': "Pressure (Pa)"},
              'velocity': {'xlabel': 'Edge', 'ylabel': "Velocity ($\mu$m/s)"},
              'pe_glc_ext': {'xlabel': 'Edge', 'ylabel': "Axial peclet number"},
              'pe_lac_ext': {'xlabel': 'Edge', 'ylabel': "Axial peclet number"}
              }
    props = {'node': ['pressure'], 'edge': ["velocity"]}

    # xaxis data:
    G = get_eucledian_lengths_wrt_origin(fs=RESULTS_SIMGRAPH_FILE)

    if measure == "pressure":
        d = nx.get_node_attributes(G, 'ed_from_source')  # eucledian distance from source node
    elif measure == "velocity":
        d = nx.get_edge_attributes(G, 'ed_from_source')

    # yaxis data
    if single:
        scan = get_static_data(fs=fpath_s, fc=fpath_c, interpolate=False)
    elif pscan:
        scan = get_pscan_static_data(f=fpath_s, interpolate=False)

    y = {}
    for key in scan.keys():
        data = scan[key]
        node_attr, edge_attr = data['node_attr'], data['edge_attr']

        if measure in props['node']:
            m = OrderedDict(zip(sorted(G.nodes), node_attr[measure]))
        elif measure in props['edge']:
            # m = OrderedDict(zip(sorted(G.edges), edge_attr[measure]))
            m = OrderedDict(zip(range(1, len(G.edges())), edge_attr[measure]))
        if sort_by_distance:
            temp = dict((v, k) for k, v in d.items())  # swaps and xdata will be position
            temp = OrderedDict(sorted(temp.items(), key=lambda kv: kv[0]))  # sorted (pos, node) / (pos, edge)
            sorted_keys = list(temp.values())
            m = OrderedDict(sorted(m.items(), key=lambda pair: sorted_keys.index(pair[0])))  # sorted(nodes, pressure)

        x = m.keys()  # xdata
        y[key] = m.values()  # ydata

    return {'x': x, 'y': y, 'xlabel':  labels[measure]['xlabel'], 'ylabel':  labels[measure]['ylabel']}


def plot_static_data(scan: List):
    """"
    plots comsol vs simgraph data
    :param data: x and multiple y values
    """
    if single:
        marker = {'simgraph': 's', 'comsol': '^'}
        color = {'simgraph': 'k', 'comsol': 'grey'}
        linewidth = {'simgraph': 0, 'comsol': 0.5}
    elif pscan:
        symbols = ['p', '^', 's', 'o', 'v', 'd', '+', 'x', "_", "1"]
        colors = ['grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey']
        widths = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        marker = OrderedDict(zip(pressure_scan, symbols))
        color = OrderedDict(zip(pressure_scan, colors))
        linewidth = OrderedDict(zip(pressure_scan, widths))

    # fig = plt.figure()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    axes = (ax1, ax2)
    for ax, data in zip(axes, scan):
        x = [str(i) for i in data['x']]  # data['x']
        y = data['y']

        for key in data['y'].keys():
            if single:
                label = key
            elif pscan:
                label = f'$\Delta$P {key} (Pa)'
            ax.plot(
                x, y[key],
                color=color[key],
                marker=marker[key],
                label=label,
                linestyle='dashed',
                linewidth=linewidth[key]
                )

        ax.set_xlabel(data['xlabel'], fontsize=12)
        ax.set_ylabel(data['ylabel'], fontsize=12)
        ax.text(
            # 0.95, 0.01, # (right, bottom)
            0.9, 0.85,
            '(' + chr(scan.index(data)+97) + ')',
            # string.ascii_uppercase[idx],
            transform=ax.transAxes,
            # weight='bold',
            size=15
        )
        xdata = [int(i) for i in x]
        # ax.set_xticks(np.arange(0, len(x)+1, 4))  #rotation=90
        # ax.set_yticks()
        leg = ax.legend(
                         loc='best',
                         ncol=1,
                         fontsize=12
                         )
        leg.get_frame().set_alpha(0.5)
        plt.setp(ax.get_xticklabels(), rotation=90)
    plt.show()

    return fig


def plot_peclet_distribution(fpath: Path, pdata: List = pressure_scan, greyscale=False):
    """
    plots peclet distributions
    reference:
    ---------
    https://stackoverflow.com/questions/43872450/matplotlib-histogram-with-multiple-legend-entries
    https://stackoverflow.com/questions/20118258/matplotlib-coloring-line-plots-by-iteration-dependent-gray-scale
    """
    scan = get_pscan_static_data(f=fpath, interpolate=False)

    fig, axes = plt.subplots(nrows=1, ncols=1)  #, figsize=(16, 5))
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    bins = 10
    if greyscale:
        colors = [(t/2.0, t/2.0, t/2.0, 0.5) for t in np.arange(0., 2., 0.4)]
        fc = dict(zip(pdata, colors))
        linestyle = dict(zip(pdata, 'solid'))

    else:
        # fc = dict(zip(pdata, [(0, 0, 1, 0.5), (1, 0, 0, 0.5), (0, 0, 0, 0.5)])) #, (1, 0, 0.5, 0.25), (0.5, 0, 0.5, 0.25)]))
        fc = dict(zip(pdata, [(0, 0, 0, 0.5), (0, 0, 0, 0.5), (0, 0, 0, 0.5)]))
        linestyle = dict(zip(pdata, ['dashed', 'solid', 'dotted']))

    for i, key in enumerate(pdata):
        data = scan[key]
        node_attr, edge_attr = data['node_attr'], data['edge_attr']
        x = edge_attr['pe_glc_ext']
        plt.hist(
            x,
            bins=bins,
            stacked=True,
            density=True,
            fc=fc[key],
            # fill=False,
            histtype='step',
            linestyle=linestyle[key],
            color='k',  # fc[key],
            linewidth=2
        )

    plt.ylabel("Density", fontsize=10)
    plt.xlabel("Axial peclet number", fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    # plt.xscale('log')

    # create legend
    # handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in fc.values()]
    handles = [Line2D([0], [0], color='k', linewidth=3, linestyle=ls) for ls in ['dashed', 'solid']]
    labels = [f"$\Delta$P = {p} (Pa)" for p in pdata]
    plt.legend(handles, labels)
    plt.show()
    return fig


def plot_conc_pscan(fpath: Path):
    """ plots conc vs nodes """

    # symbols = ['+', '^', 's', 'p', 'v', 'd', 'o', 'x', "_", "1"]
    # widths = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    # marker = OrderedDict(zip(pressure_scan, symbols))
    # linewidth = OrderedDict(zip(pressure_scan, widths))

    fig = plt.figure()
    time, scan = get_pscan_dynamic_data(f=fpath, resample=False)
    print(scan.keys())
    # nodes to plot (inlet, center, outlet) (34, 4 , 39)
    nodes = [27, 11, 38]
    color = OrderedDict(zip(nodes, ['r', 'g', 'b']))
    linestyle = OrderedDict(zip(scan.keys(), ['dashed', 'solid']))

    data = {}
    for idx, s in enumerate(scan):
        species = scan[s]["glc_ext"]
        data[s] = pd.DataFrame(species)
    for key in data.keys():
        print(key)
        for n in nodes:
            plt.plot(
                time, data[key][n],
                color=color[n],
                # marker=marker[key],
                markersize=5,
                label=key,
                linestyle=linestyle[key],
                # linewidth=linewidth[key]
            )
    plt.xlabel('Time (s)', fontsize=10)
    plt.ylabel('Concentration (mM)', fontsize=10)
    plt.show()
    return fig


def load_G_H():
    G = create_graph(attributes=True)
    input = get_graph(sheet="test11")
    H = create_graph(output=input, derived_attr=True)
    return G, input, H


class ManuscriptFigures:
    def figures(self):
        # self.draw_graph3d_vedo()
        self.plot_3dhistogram()
        # self.plot_morphological_params()
        # self.plot_qdot()
        # self.Fig1()  # static properties single run comparison of comsol vs. simgraph
        # self.Fig2() # static properties pscan run
        self.Fig3()  # pe distribution
        # self.Fig4()  # plot of pressure scan at diferent pe
        # self.Fig5()
        # self.Fig6()  # conc plot of simgraph vs comsol
        # self.Fig7()
        # self.Fig8()
        # self.Fig9()  #line plot of pressure and velocity
        # self.Fig10()  # methods figure 2
        # self.Fig11()  # volume interpolation plot for glc, lac
        # self.Fig12()  # probe line plot for glc and lac exchange
        # self.Fig13()  # plots steady-state/rise times of simgraph vs. comsol
        # self.Fig14()  # plots steady-state/rise times of tumor design 1 vs. design 2
        # self.Fig15()    # plots bv cell glucose and lactate time-varying concentration in 3D vasculature
        # self.Fig16()    # glucose uptake and latate release flux
        # self.Fig17()
        # self.Fig18()
        # self.Fig19()

    def plot_3dhistogram(self):
        """
        :param data: Attribute data to plot
        :param label: xaxis label
        :param binwidth:
        :param density:
        :return:
        reference:
        https://stackoverflow.com/questions/6871201/plot-two-histograms-on-single-chart-with-matplotlib

        example:
        import random
        import numpy
        from matplotlib import pyplot

        x = [random.gauss(3, 1) for _ in range(400)]
        y = [random.gauss(4, 2) for _ in range(400)]

        bins = numpy.linspace(-10, 10, 100)

        pyplot.hist(x, bins, alpha=0.5, label='x')
        pyplot.hist(y, bins, alpha=0.5, label='y')
        pyplot.legend(loc='upper right')
        pyplot.show()
        """
        fc = dict(zip(test_cases, ['b', 'r', 'k', 'g']))
        linestyle = dict(zip(test_cases, ['solid', (0, (5, 10)), 'dotted', 'solid']))

        fig, axes = plt.subplots(nrows=2, ncols=1)  # , figsize=(16, 5))
        # fig.subplots_adjust(wspace=0.3, hspace=0.3)

        data = {}
        for sheet in test_cases:
            data[sheet] = get_graph(sheet=sheet)

        param = {'d': {'bw': 5, 'label': 'Diameter', 'figlabel': 'e'}, 'l': {'bw': 10, 'label': 'Length', 'figlabel': 'f'}}
        for i, p in enumerate(param):
            for idx, key in enumerate(data):
                x = data[key][p]
                ax = axes[i]
                ax.hist(
                    x,
                    # bins=10,
                    stacked=False,
                    density=True,
                    histtype='step',
                    linestyle=linestyle[key],
                    color=fc[key],
                    linewidth=2
                )
            # create legend
            handles = [Line2D([0], [0], color=c, linewidth=3, linestyle=ls) for c, ls in zip(fc.values(), linestyle.values())]
            labels = ["islet", "tumor design 1", "tumor design 2", "mesentery"]
            ax.legend(handles, labels, prop={"size":10})
            ax.set_ylabel("Density", fontsize=12)
            ax.set_xlabel(f"{param[p]['label']} ($\mu$m)", fontsize=12)
            ax.set_xscale('log')
            ax.text(
                0.008, top,
                # '(' + chr(i + 97) + ')',
                '(' + param[p]['figlabel'] + ')',
                transform=ax.transAxes,
                size=12
            )
        plt.show()
        f = os.path.join(RESULTS_PLOTS_DIR, "morphology.svg")
        fig.savefig(f, transparent=True, dpi=600, bbox_inches="tight")

    def draw_graph3d_vedo(self):
        """
        :return:
        """
        # model , label : font size

        font_size = {"test9_default": 35, "test10_default": 35, "test11_default": 5, "test24_default": 125}
        settings.windowSplittingPosition = 0.4
        vp = Plotter(shape='3/1', interactive=True, sharecam=False)

        for i, sheet in enumerate(test_cases):
            input = get_graph(sheet=sheet)
            source = int(input['source'])
            target = int(input['target'])
            t = {source: 'In', target: 'Out'}  # terminals inlet exit
            c = {source: (0, 0, 0), target: (0, 0, 0)}

            G = create_graph(output=input)
            nodes = sorted(G.nodes())

            # create labels
            l = []
            color = []
            for n in nodes:
                if n in t.keys():
                    l.append(t[n])
                    color.append(c[n])
                else:
                    l.append('')
                    color.append((128, 128, 128))

            # create vtk polydata
            pos = dict(zip(input['nodes'], input['xyz']))
            lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
            pts = Points(input['xyz'], r=8, c=color).rotateX(0)
            edg = Lines(lines).lw(5).rotateX(0).c("gray")

            scale = font_size[sheet]
            labs = pts.labels(l, scale=scale, font='VictorMono', c='k').addPos(0.1, 1, 0.5) #FIXME: revert to (0.05,0,0.5)
            f = os.path.join(RESULTS_PLOTS_DIR, f'{sheet}.png')

            if True:
                plt = Plotter(N=1)  # , size=(800, 600)
                if sheet == "test11_default":
                    axes = Axes(
                        edg, xtitle='y (\mum)', ytitle='z (\mum)', ztitle='x (\mum)', xyGrid=True, zxGrid=True, yzGrid=True,
                        xLabelSize=0.022, yLabelSize=0.022, zLabelSize=0.022, xTitleSize=0.03, yTitleSize=0.03,
                        zTitleSize=0.03
                    )
                elif sheet == "test9_default":
                    axes = Axes(
                        edg, xyGrid=True, zxGrid=True, yzGrid=True, xtitle='x (\mum)', ytitle='y (\mum)', ztitle='z (\mum)',
                        xLabelSize=0.022, yLabelSize=0.022, zLabelSize=0.022, xTitleSize=0.03, yTitleSize=0.03,
                        zTitleSize=0.03, xTitleOffset=-0.01, yTitleOffset=-0.05,
                        xTitlePosition=1.05, yTitlePosition=1.05,
                        )
                elif sheet == "test10_default":
                    axes = Axes(
                        edg, xyGrid=True, zxGrid=True, yzGrid=True, xtitle='x (\mum)', ytitle='y (\mum)',
                        ztitle='z (\mum)',
                        xLabelSize=0.022, yLabelSize=0.022, zLabelSize=0.022, xTitleSize=0.03, yTitleSize=0.03,
                        zTitleSize=0.03, xTitleOffset=-0.01, yTitleOffset=-0.05,
                        xTitlePosition=1.05, yTitlePosition=1.1,
                    )
                elif sheet == "test24_default":
                    axes = Axes(
                        edg,
                        xyGrid=True, zxGrid=True, yzGrid=True, xTitlePosition=0.15, yTitlePosition=1.02,
                        xtitle='x (\mum)', ytitle='y (\mum)', ztitle='z (\mum)',
                        yLabelRotation=-90, xLabelRotation=-90, xTitleRotation=0, yTitleRotation=-180, xTitleOffset=-0.06, yTitleOffset=0.001,
                        xTitleSize=0.022, yTitleSize=0.022, zTitleSize=0.022
                    )
                plt.show(
                    pts,
                    edg,
                    labs,
                    axes,
                    bg2='w',
                    interactive=True,
                    # zoom=1.2,
                    at=0,
                    # roll=90
                ).screenshot(f)
                interactive()
            if False:
                axes = Axes(edg, xyGrid=True, zxGrid=True, yzGrid=True)  # zLabelRotation=-45)
                vp.show(
                    pts,
                    edg,
                    labs,
                    axes,
                    bg2='w',
                    at=i
                )

    def plot_morphological_params(self):
        """
        generates histogram plots of length and diameter distribution
        :return:

        Reference:
        ----------
        enumerate subplot: https://stackoverflow.com/questions/22508590/enumerate-plots-in-matplotlib-figure
        enumerate subplot: https://stackoverflow.com/a/25544329/8281509
        position text label: https://stackoverflow.com/questions/8482588/putting-text-in-top-left-corner-of-matplotlib-plot
        """

        data = {}
        for sheet in test_cases:
            data[sheet] = get_graph(sheet=sheet)

        fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(nrows=2, ncols=4, figsize=(16, 12))
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        axes = (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)
        param = {'d': {'bw': 5, 'label': 'Diameter'}, 'l': {'bw': 10, 'label': 'Length'}}

        k = 0

        for p in param.keys():
            for idx, test_case in enumerate(data):
                bw = param[p]['bw']
                label = param[p]['label']
                ax = axes[k]
                if idx == 0:
                    ax.set_ylabel('Frequency', fontsize=12)

                ax.set_xlabel(f"{label} ($\mu$m)", fontsize=12)
                ax.text(
                    right, top,
                    '('+fig_labels[test_case]+')',
                    # string.ascii_uppercase[idx],
                    transform=ax.transAxes,
                    # weight='bold',
                    size=20
                )
                plot_histogram(ax, data=data[test_case][p], binwidth=bw)
                k = k + 1
        plt.show()
        f = os.path.join(RESULTS_PLOTS_DIR, "morphology.svg")
        fig.savefig(f, transparent=True, dpi=600, bbox_inches="tight")

    def plot_qdot(self):
        weights = []
        for test_case in test_cases:
            fpath = os.path.join(RESULTS_SIMGRAPH_DIR, f"{test_case}_bc1_v1_c0.mat")
            scan = get_static_data(fs=fpath, interpolate=False)
            for idx, key in enumerate(scan):
                data = scan[key]
                node_attr, edge_attr = data['node_attr'], data['edge_attr']
                weights.append(edge_attr['qdot'])

        fig = draw_weighted_connectivity_matrix(weights=weights, test_cases=test_cases)

        if pb:
            f = "qdot_pb"
        else:
            f = "qdot"

        fpath = os.path.join(RESULTS_PLOTS_DIR, f'{f}.svg')
        fig.savefig(fpath, transparent=True, bbox_inches='tight', dpi=600)
        plt.show()

    def Fig1(self):
        """ static properties for single run """
        for test_case in test_cases:
            fs = os.path.join(RESULTS_SIMGRAPH_DIR, f"{test_case}_bc1_v1_c0.mat")
            fc = os.path.join(RESULTS_COMSOL_DIR, f"{test_case}_bc1_v1_c0.mat")
            print(fs)
            # pressure
            pressure_data = get_plot_data(
                fpath_s=fs,
                fpath_c=fc,
                measure="pressure",
                sort_by_distance=False
            )

            # velocity
            velocity_data = get_plot_data(
                fpath_s=fs,
                fpath_c=fc,
                measure="velocity",
                sort_by_distance=False
            )

            fig = plot_static_data(scan=[pressure_data, velocity_data])
            f = os.path.join(RESULTS_PLOTS_DIR, f"{test_case}_pressure_velocity_scatter.svg")
            fig.savefig(f, transparent=False, dpi=600, bbox_inches="tight")

    def Fig2(self):
        """ static properties for pscan run """
        print(RESULTS_PSCAN_FILE)
        pressure_data = get_plot_data(
            fpath_s=RESULTS_PSCAN_FILE,
            fpath_c=None,
            measure="pressure",
            sort_by_distance=False
        )

        # velocity
        velocity_data = get_plot_data(
            fpath_s=RESULTS_PSCAN_FILE,
            fpath_c=None,
            measure="velocity",
            sort_by_distance=False
        )

        fig = plot_static_data(scan=[pressure_data, velocity_data])
        f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_pressure_velocity_scatter_pscan.svg")
        fig.savefig(f, transparent=True, dpi=600, bbox_inches="tight")

    def Fig3(self):
        """ plot 2D histogram of axial peclet numbers for species = glc_ext """
        fig = plot_peclet_distribution(
            fpath=RESULTS_PSCAN_FILE,
            pdata=[20, 200]
        )
        f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_pe_pscan.svg")
        fig.savefig(
            f,
            transparent=True,
            dpi=600,
            bbox_inches="tight"
        )

    def Fig4(self):
        """ plot conc of pscan at different peclet numbers and different nodal ppositions of the vasculature"""

        # symbols = ['+', '^', 's', 'p', 'v', 'd', 'o', 'x', "_", "1"]
        # widths = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        # marker = OrderedDict(zip(pressure_scan, symbols))
        # linewidth = OrderedDict(zip(pressure_scan, widths))
        # gimg = os.path.join(RESULTS_PLOTS_DIR, f'{test}.png')

        # ax = plt.subplot(111)
        fig, ax = plt.subplots()
        time, scan = get_pscan_dynamic_data(f=RESULTS_PSCAN_FILE, resample=False)
        # nodes to plot (inlet, center, outlet) (34, 4 , 39)
        nodes = [27, 5, 38] #5>11
        color = OrderedDict(zip(nodes, ['r', 'g', 'b']))
        print(scan.keys())
        linestyle = OrderedDict(zip(pressure_scan, ['dashed', 'solid']))

        # --------------------------------------------------------------------------------------------------------------
        # main plot
        # --------------------------------------------------------------------------------------------------------------

        data = {}
        for idx, s in enumerate(scan):
            species = scan[s]["glc_ext"]
            data[s] = pd.DataFrame(species)
        for key in data.keys():
            if key in pressure_scan:
                for n in nodes:
                    ax.plot(
                        time, data[key][n],
                        color=color[n],
                        # marker='*',#marker[key],
                        # markersize=5,
                        label=key,
                        linestyle=linestyle[key],
                        # linewidth=linewidth[key]
                    )
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_ylabel('Concentration (mM)', fontsize=12)

        # add inset
        # arr_img = plt.imread(gimg)
        # im = OffsetImage(arr_img, zoom=72. / 120)
        # ab = AnnotationBbox(im, (0.75, 0.25), xycoords='axes fraction', frameon=False)
        # ax.add_artist(ab)

        # --------------------------------------------------------------------------------------------------------------
        # inset plot
        # --------------------------------------------------------------------------------------------------------------

        rect = [.45, 0.24, .5, .4]
        bbox = Bbox.from_bounds(*rect)
        inax = fig.add_axes(bbox)
        inax.axis('on')
        inax.set_facecolor('none')
        # inax.grid('on', color='k')

        scan = get_pscan_static_data(f=RESULTS_PSCAN_FILE, interpolate=False)
        bins = 10
        fc = dict(zip(pressure_scan, [(0, 0, 0, 0.5), (0, 0, 0, 0.5), (0, 0, 0, 0.5)]))
        linestyle = dict(zip(pressure_scan, ['dashed', 'solid', 'dotted']))

        for i, key in enumerate(pressure_scan):
            data = scan[key]
            node_attr, edge_attr = data['node_attr'], data['edge_attr']
            x = edge_attr['pe_glc_ext']
            inax.hist(
                x,
                bins=bins,
                stacked=True,
                density=True,
                fc=fc[key],
                # fill=False,
                histtype='step',
                linestyle=linestyle[key],
                color='k',  # fc[key],
                linewidth=2
            )

        inax.set_ylabel("Density", fontsize=10)
        inax.set_xlabel("Axial peclet number", fontsize=10)

        fig.show()
        f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_conc_pscan.svg")
        ax.figure.savefig(
            f,
            transparent=True,
            dpi=600,
            bbox_inches="tight"
        )
        # f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_conc_pscan.svg")
        # ax.figure.savefig(f, transparent=True, dpi=600, bbox_inches="tight")
        # plt.show()

    def Fig5(self):
        """plot fluxes of pscan"""
        fig = plt.figure()
        data = io_utils_py._get_file(RESULTS_PSCAN_FILE, spreadsheet=False)
        glcim, lacex, delp = [], [], []
        for s in data['variable']:
            glcim.append(s.glcim_net)
            lacex.append(s.lacex_net)
            delp.append(s.delP)

        plt.plot(
            delp,
            glcim,
            color='k'
        )
        plt.xlabel('Time (s)')
        plt.ylabel('Concentration (mM)')
        plt.show()
        f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_flux_pscan.svg")
        fig.savefig(f, transparent=True, dpi=600, bbox_inches="tight")

    def Fig6(self):
        """
        plot conc for comparing simgraph vs. comsol
        inset figure is loaded from a file
        references:
        ----------
        https://stackoverflow.com/questions/66540026/insert-an-svg-image-in-matplotlib-figure
        https://stackoverflow.com/questions/12161559/plot-svg-within-matplotlib-figure-embedded-in-wxpython
        """
        gimg = os.path.join(RESULTS_PLOTS_DIR, f'{test}_inset.png')
        fs = os.path.join(RESULTS_SIMGRAPH_DIR, f"{test}_bc1_v1_c0.mat")
        fc = os.path.join(RESULTS_COMSOL_DIR, f"{test}_bc1_v1_c0.mat")
        time, scan = get_dynamic_data(fs=fs, fc=fc)

        # --------------------------------------------------------------------------------------------------------------
        # settings and labels
        # --------------------------------------------------------------------------------------------------------------
        linestyle = {'simgraph': 'solid', 'comsol': 'dashed'}
        if test == "test11_default":
            nodes = [27, 5, 38]
            colors = ['r', 'g', 'b']
        elif test == "test24_default":
            # nodes = [79, 19, 531, 13, 606, 577, 151, 311]  # neighbors of 13 (53,606,577)
            # colors = ['r', 'c', 'y', 'k', 'brown', 'g', 'm', 'b']
            nodes = [79, 19, 13, 151, 311]  # highlight nodes
            colors = ['r', 'c', 'g', 'm', 'b']
        color = OrderedDict(zip(nodes, colors))

        # --------------------------------------------------------------------------------------------------------------
        # main plot
        # --------------------------------------------------------------------------------------------------------------
        ax = plt.subplot(111)
        for key in scan.keys():
            species = pd.DataFrame(scan[key]["glc_ext"])
            print(species)
            for n in nodes:
                ax.plot(
                    time, species[n-1],
                    color=color[n],
                    markersize=5,
                    label=key,
                    linestyle=linestyle[key],
                )
        ax.set_xlabel('Time (s)', fontsize=10)
        ax.set_ylabel('Concentration (mM)', fontsize=10)

        # --------------------------------------------------------------------------------------------------------------
        # inset plot
        # --------------------------------------------------------------------------------------------------------------
        if True:
            draw_graph3d_vedo(
                point_r=10,
                lw=2,
                points=False,
                sphere=True,
                highlight=True,
                line_alpha=1,
                fout=f"{test}_inset.png",
                default_axes=False,
                interact=True,
                sphere_r=30, #test24
                label_scale=100, #test24
                label_xoffset=80 #test24
            )

        # G = get_eucledian_lengths_wrt_origin(fs=RESULTS_SIMGRAPH_FILE, origin=34)
        # d = nx.get_node_attributes(G, "ed_from_source")
        if False:
            arr_img = plt.imread(gimg)
            im = OffsetImage(arr_img, zoom=0.45)
            ab = AnnotationBbox(im, (1, 0), xycoords='axes fraction', box_alignment=(1.1, -0.1), frameon=False)
            ax.add_artist(ab)
        plt.show()

        f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_conc_compare.svg")
        ax.figure.savefig(
            f,
            transparent=True,
            dpi=600,
            bbox_inches="tight"
        )

    def Fig7(self):
        """ plot conc fo pscan"""
        fs = os.path.join(RESULTS_SIMGRAPH_DIR, f"{test}_pb_bc1_v1_c0.mat")
        fc = os.path.join(RESULTS_COMSOL_DIR, f"{test}_pb_bc1_v1_c0.mat")
        time, scan = get_dynamic_data(fs=fs, fc=fc)

        # --------------------------------------------------------------------------------------------------------------
        # settings and labels
        # --------------------------------------------------------------------------------------------------------------
        linestyle = {'simgraph': 'solid', 'comsol': 'dashed'}
        nodes = [27, 11, 38]
        colors = ['r', 'g', 'b']
        color = OrderedDict(zip(nodes, colors))

        # --------------------------------------------------------------------------------------------------------------
        # inset plot
        # --------------------------------------------------------------------------------------------------------------

        vplt = draw_graph3d_vedo(
            point_r=10,
            lw=2,
            points=False,
            sphere=True,
            highlight=True,
            line_alpha=1,
            rotatex=0,
        )
        np_pic = Picture(screenshot(returnNumpy=True, scale=1))
        vplt.close()

        # --------------------------------------------------------------------------------------------------------------
        # main plot
        # --------------------------------------------------------------------------------------------------------------
        ax = plt.subplot(111)
        for key in scan.keys():
            species = pd.DataFrame(scan[key]["glc_ext"])
            for n in nodes:
                ax.plot(
                    time, species[n],
                    color=color[n],
                    markersize=5,
                    label=key,
                    linestyle=linestyle[key],
                )
        ax.set_xlabel('Time (s)', fontsize=10)
        ax.set_ylabel('Concentration (mM)', fontsize=10)

        f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_conc_compare.svg")
        ax.figure.savefig(
            f,
            transparent=True,
            dpi=600,
            bbox_inches="tight"
        )

        # ------------------------------------------------------------
        # Build plot: exponential
        # ------------------------------------------------------------
        x = np.arange(0, 4, 0.1)
        y = 3 * np.exp(-x)

        plt1 = plot(
            x, y,
            xtitle='time in seconds',
            ytitle='some function [a.u.]',
        )

        np_pic.scale(0.025).x(2).y(1.)

        show(plt1, np_pic, size=(800, 600), zoom=1.2)

    # ------------------------------------------------------------------------------------------------------------------
    # figure8
    # ------------------------------------------------------------------------------------------------------------------
    def Fig8(self):
        """ plot conc fo pscan using 3d plot generated in matplotlib
        reference:
        """
        fs = os.path.join(RESULTS_SIMGRAPH_DIR, f"{test}_bc1_v1_c0.mat")  # FIXME:revert to _pb for test11
        fc = os.path.join(RESULTS_COMSOL_DIR, f"{test}_bc1_v1_c0.mat")
        time, scan = get_dynamic_data(fs=fs, fc=fc)
        print(scan.keys())
        # --------------------------------------------------------------------------------------------------------------
        # settings and labels
        # --------------------------------------------------------------------------------------------------------------
        G = create_graph()
        terminals = [335, 352]  # FIXME: revert to [34, 39] for test11 # source & target
        hnodes = [27, 11, 38]  # FIXME: revert to [40, 89, 123, 161, 312, 352]  for test24 # highlight nodes
        terminal_labels = dict(zip(terminals, ['In', 'Out']))
        hnode_color = OrderedDict(zip(hnodes + terminals, ['r', 'g', 'b', 'y', 'c', 'm', 'k', 'k']))  # FIXME: revert to ['r', 'g', 'b', 'k', 'k'] for test11
        hnode_alpha = dict.fromkeys(hnodes + terminals, 1)
        label = dict.fromkeys(G.nodes(), '')
        color = dict.fromkeys(G.nodes(), "grey")
        alpha = dict.fromkeys(G.nodes(), 0.3)
        label.update(terminal_labels)
        color.update(hnode_color)
        alpha.update(hnode_alpha)
        linestyle = {'simgraph': 'solid', 'comsol': 'dashed'}

        # --------------------------------------------------------------------------------------------------------------
        # main plot
        # --------------------------------------------------------------------------------------------------------------
        fig, ax = plt.subplots()
        for key in scan.keys():
            species = pd.DataFrame(scan[key]["glc_ext"])
            for n in hnodes:
                ax.plot(
                    time, species[n],
                    color=color[n],
                    marker='*',
                    markersize=5,
                    label=key,
                    linestyle=linestyle[key],
                )
        ax.set_xlabel('Time (s)', fontsize=10)
        ax.set_ylabel('Concentration (mM)', fontsize=10)

        # --------------------------------------------------------------------------------------------------------------
        # inset plot
        # --------------------------------------------------------------------------------------------------------------

        # rect = [.3, 0.1, .8, .7]
        rect = [.5, 0.24, .4, .3]
        bbox = Bbox.from_bounds(*rect)
        inax = fig.add_axes(bbox) #, projection='3d')
        inax.axis('on')
        # inax.view_init(azim=90) #, elev=0)
        inax.set_facecolor('none')
        # inax.grid('on', color='k')

        seen = set()
        for i, j in G.edges():
            c1 = nx.get_node_attributes(G, 'pos')[i]
            c2 = nx.get_node_attributes(G, 'pos')[j]
            x = np.stack((c1, c2))
            inax.plot(*x.T, color='grey')

            if i not in seen:
                inax.scatter(*x[0], color=color[i], alpha=alpha[i], s=15)
                inax.text(*x[0], label[i], size=10, color='k')
                seen.add(i)
            if j not in seen:
                inax.scatter(*x[1], color=color[j], alpha=alpha[j], s=15)
                inax.text(*x[1], label[j], size=10, color='k')
                seen.add(j)

            inax.set_xlabel('x')
            inax.set_ylabel('y')
            inax.set_zlabel('z')
            inax.set_xticklabels([])
            inax.set_yticklabels([])
            inax.set_zticklabels([])
            inax.tick_params(axis="x", colors="white")
            inax.tick_params(axis="y", colors="white")
            inax.tick_params(axis="z", colors="white")
            inax.set_xlabel('x', labelpad=-10, loc='right')
            inax.set_ylabel('y', labelpad=-10, loc='top')
            inax.set_zlabel('z', labelpad=-10)
            inax.w_xaxis._axinfo.update({'grid': {'color': 'lightgrey', 'linewidth': 0.8, 'linestyle': '-'}})
            inax.w_yaxis._axinfo.update({'grid': {'color': 'lightgrey', 'linewidth': 0.8, 'linestyle': '-'}})
            inax.w_zaxis._axinfo.update({'grid': {'color': 'lightgrey', 'linewidth': 0.8, 'linestyle': '-'}})


        fig.show()
        f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_conc_compare.svg")
        ax.figure.savefig(
            f,
            transparent=True,
            dpi=600,
            bbox_inches="tight"
        )

    # ------------------------------------------------------------------------------------------------------------------
    # figure9
    # ------------------------------------------------------------------------------------------------------------------
    def Fig9(self):
        """
        compare velocity simgraph for a pressure scan
        :param cname: colormap name
        :param fpath_c:
        :param fpath_s:
        :param pts:
        :param edg:
        :return:
        """
        def plot_pressure():
            plt1 = Plotter(interactive=False)
            vmin, vmax = get_bounds_static(scan=scan, param='pressure', attr='node_attr')

            # manually build a lookup table for mapping colors
            lut = buildLUT([
                (vmin, 'dg',),
                (vmax, 'red'),
            ],
                vmin=vmin, vmax=vmax,
                interpolate=True,
                belowColor='yellow',
            )
            pts1 = pts.clone().pointSize(5).pointColors(node_attr['pressure'], cmap='bwr', vmin=vmin, vmax=vmax)
            pts1.addScalarBar3D(
                title='(a) pressure (Pa)',
                titleFont='VictorMono',
                labelFont='VictorMono',
                c='k',
            )
            pts1.scalarbar.addPos(1, 0, 0).useBounds()
            edg1 = edg.clone().pointColors(edge_attr['pressure'], cmap='bwr')
            for a in pts1.scalarbar.unpack():
                if a and a.name == 'Text': a.scale(1.5)

            axes = Axes(edg1, xtitle='y', ytitle='z', ztitle='x', xyGrid=True, zxGrid=True, yzGrid=True)  # zLabelRotation=-45)
            # plt1.add(pts1, edg1, axes)
            f = os.path.join(RESULTS_PLOTS_DIR, f'{test}_pressure.png')
            pltp = plt1.show(pts1, edg1, axes).screenshot(f)
            return pltp

        def plot_velocity():
            """ line plot of velocity: 3D visualization"""
            vmin, vmax = get_bounds_static(scan=scan, param='velocity', attr='edge_attr', last=False)

            plt2 = Plotter(interactive=False)
            # manually build a lookup table for mapping colors
            lut = buildLUT([
                (vmin, 'b'),
                ((vmin+vmax)/2, 'w'),
                (vmax+2, 'r'),
            ],
                vmin=vmin, vmax=vmax+2,
                interpolate=True,
                aboveColor='dr',
            )
            color = []
            for v in edge_attr['velocity']:
                rgb = [0, 0, 0]
                lut.GetColor(v, rgb)
                color.append(getColorName(rgb))
            edge_w = OrderedDict(zip(sorted(G.edges), color))
            nx.set_edge_attributes(G, edge_w, 'color')

            nodes = pts.points()
            arrsv = []
            for i, j, attr in G.edges(data=True):
                dv = nodes[j - 1] - nodes[i - 1]
                cn = (nodes[j - 1] + nodes[i - 1]) / 2
                v = dv / mag(dv) * 3.5
                ar = Arrow(cn - v, cn + v, s=.1, c=G[i][j]["color"])  # attr['color']
                arrsv.append(ar)
            edg2 = edg.clone().cmap(lut, edge_attr['velocity'], on='cells').c('grey')
            edg2.addScalarBar3D(
                title='(b) velocity (\mum/s)',
                titleFont='VictorMono',
                labelFont='VictorMono',
                c='k',
            )
            edg2.scalarbar.addPos(0.1, 0, 0).useBounds()
            for a in edg2.scalarbar.unpack():
                if a and a.name == 'Text': a.scale(1.4)

            axes = Axes(edg2, xtitle='y', ytitle='z', ztitle='x', xyGrid=True, zxGrid=True, yzGrid=True)  # zLabelRotation=-45)
            # plt2.add(pts, edg2, arrsv, axes)
            f = os.path.join(RESULTS_PLOTS_DIR, f'{test}_velocity.png')
            pts2 = pts.clone().pointSize(5).c('gray')
            pltv = plt2.show(pts2, edg2, arrsv, axes).screenshot(f)

            return pltv


        input = get_graph(sheet=None)
        # create vtk polydata
        pos = dict(zip(input['nodes'], input['xyz']))
        lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]

        # define points and lines for vtk
        pts = Points(input['xyz'], r=10).lighting('off').rotateX(-90)
        edg = Lines(lines).lw(5).rotateX(-90)

        scan = get_static_data(fs=RESULTS_SIMGRAPH_FILE, fc=None)
        data = scan["simgraph"]
        node_attr, edge_attr = data['node_attr'], data['edge_attr']

        pltp = plot_pressure()
        pltv = plot_velocity()
        f = os.path.join(RESULTS_PLOTS_DIR, f'{test}_p_v.png')
        plt12 = Plotter(N=2, interactive=True, sharecam=True)
        plt12.show(pltp.actors, at=0, zoom=1.2)
        plt12.show(pltv.actors, at=1, zoom=1.2).screenshot(f)

    def Fig10(self):
        """ creates methods figure displaying the cylindrical and spherical volume elements"""
        if False:
            G, input, H = load_G_H()
            G_nodes = list(G.nodes())
            source = int(input['source'])
            target = int(input['target'])

            # node shapes and attributes
            s_nodes = [n for n in G_nodes if n not in [source, target]] #sphere
            c_nodes = list(set(H.nodes) - set(s_nodes))
            pos = nx.get_node_attributes(H, 'pos')
            node_r = nx.get_node_attributes(H, 'node_r')
            node_l = nx.get_node_attributes(H, 'node_l')
            # create vtk object
            nx_lines = [(pos[t], pos[h]) for t, h in G.edges()]
            edg = Lines(nx_lines).lw(5)

            pts = []
            for n in H.nodes():
                if n in s_nodes:
                    pts.append(Sphere(r=node_r[n]*0.5, alpha=0.8).pos(pos[n]))
                if n in c_nodes:
                    cyl1 = Cylinder(r=node_r[n]*0.3, height=node_l[n]*.5, cap=False, alpha=1, c='grey').pos(pos[n])
                    cyl2 = Cylinder(r=node_r[n]*0.5, height=node_l[n]*0.5, cap=False, alpha=0.2, c='grey').pos(pos[n])
                    cyl = merge(cyl1, cyl2)
                    pts.append(cyl)
            show(pts, edg, axes=True, interactive=True, bg='w', title='plot').screenshot("vol")

        G = nx.OrderedGraph()
        G.add_edges_from(ebunch_to_add=[(0, 12), (12, 1), (1, 2), (2, 3), (2, 4), (3, 10), (10, 5), (4, 11), (11, 6)])
        # nx_pos = {1: [0, 0, 0], 2: [12.4, 0, 0], 3: [23.14, 6.2, 0], 4: [23.14, -6.2, 0]} #30 degree elevation
        # nx_pos = {1: [0, 0, 0], 2: [12.4, 0, 0], 3: [21.17, 8.77, 0], 4: [21.17, -8.77, 0]} #45 degree elevation
        nx_pos = {0: [-12.4, 0, 0], 1: [0, 0, 0], 2: [12.4, 0, 0], 3: [18.60, 10.74, 0], 4: [18.60, -10.74, 0], 5: [24.8, 21.47, 0], 6: [24.8, -21.47, 0],
                  7: [9.3, 0, 0], 8: [13.95, 2.69, 0], 9: [13.95, -2.69, 0], 10: [21.7, 16.10, 0], 11: [21.7, -16.1, 0],
                  12: [-6.2, 0, 0]} #60 degree elevation
        # pos = {1: [[-6.2, 0, 0], [6.2, 0, 0]], 2: [12.4, 0, 0], 3: [[17.77, 3.1, 0], [28.5, 9.3, 0]], 4: [[17.77, -3.1, 0], [28.5, -9.3, 0]]} #30 degree elevation
        # pos = {1: [[-6.2, 0, 0], [6.2, 0, 0]], 2: [12.4, 0, 0], 3: [[16.78, 4.38, 0], [25.9, 13.15, 0]], 4: [[16.78, -4.38, 0], [25.9, -13.15, 0]]} #45 degree elevation
        pos1 = {1: [[-6.2, 0, 0], [6.2, 0, 0]], 2: [12.4, 0, 0], 3: [[15.5, 5.37, 0], [21.7, 16.10, 0]], 4: [[15.5, -5.37, 0], [21.7, -16.1, 0]],
                0: [[-18.6, 0, 0], [-6.2, 0, 0]], 5: [[21.7, 16.10, 0], [27.9, 26.85, 0]], 6: [[21.7, -16.10, 0], [27.9, -26.85, 0]]} #60 degree elevation
        r1 = {0: 2.6, 1: 2.6, 2: 4.55, 3: 2.45, 4: 2.6, 5: 2.45, 6: 2.6, 7: 2.6, 8: 2.45, 9: 2.6}  #FIXME: revert 2: 4.55
        r2 = {1: 8.8, 3: 8.65, 4: 8.8}
        color = {1: 'lightsteelblue', 2: 'r', 3: 'darkseagreen', 4: 'grey', 7: 'lightsteelblue', 8: 'darkseagreen',
                 9: 'grey'}
        nx.set_node_attributes(G, pos1, 'pos')
        nx.set_node_attributes(G, r1, 'r1')
        nx.set_node_attributes(G, r2, 'r2')

        plt = Plotter(shape=(1, 4), interactive=True, sharecam=0)
        # --------------------------------------------------------------------------------------------------------------
        # subfigure 0
        # --------------------------------------------------------------------------------------------------------------
        pos2 = {1: [[-12.4, 0, 0], [12.4, 0, 0]], 3: [[12.4, 0, 0], [24.8, 21.47, 0]], 4: [[12.4, 0, 0], [24.8, -21.47, 0]]}
        nx_lines = [(nx_pos[t], nx_pos[h]) for t, h in G.edges()]
        edg2 = Lines(nx_lines).lw(2)

        pts2 = []
        for n in G.nodes():
            if n in [0, 2, 5, 6]:
                pts2.append(Sphere(r=0.5, pos=nx_pos[n], c='k'))
            elif n in [1, 3, 4]:
                pts2.append(Cylinder(r=r1[n], height=12.4, cap=True, pos=pos2[n], c=color[n]).alpha(0.2))

        plt.show(pts2, edg2, axes=False, interactive=True, bg='w', at=0)

        # --------------------------------------------------------------------------------------------------------------
        # subfigure 1
        # --------------------------------------------------------------------------------------------------------------
        pos3 = {7: [[6.2, 0, 0], [12.4, 0, 0]], 8: [[12.4, 0, 0], [15.5, 5.37, 0]], 9: [[12.4, 0, 0], [15.5, -5.37, 0]]}
        pts3, edg3 = [], []
        for i, j in G.edges():
            p1 = nx_pos[i]
            p2 = nx_pos[j]
            if (i, j) in [(0, 12), (10, 5), (11, 6)]:
                edg3.append(DashedLine([p1, p2], spacing=0.3).lw(2))
            else:
                edg3.append(Line([p1, p2]).lw(2))

        for n in range(0, 10):
            if n in [7, 8, 9]:
                pts3.append(Sphere(r=0).pos(nx_pos[n]).alpha(0.4))
                pts3.append(Cylinder(r=r1[n], height=6.2, cap=True, pos=pos3[n], c=color[n]).alpha(0.4))
            elif n in [0, 1, 3, 4, 5, 6]:
                pts3.append(Sphere(r=0.5, c='k').pos(nx_pos[n]).alpha(1))
                if n in [1, 3, 4]:
                    pts3.append(Cylinder(r=r1[n], height=12.4, cap=True, pos=pos1[n], c=color[n]).alpha(0.2))
            elif n == 2:
                pts3.append(Sphere(r=0.5, c='k').pos(nx_pos[n]).alpha(1))

        plt.show(pts3, edg3, axes=False, interactive=True, bg='w', at=1)

        # --------------------------------------------------------------------------------------------------------------
        # subfigure 2
        # --------------------------------------------------------------------------------------------------------------
        pts4, edg4 = [], []
        sph1 = []
        for i, j in G.edges():
            p1 = nx_pos[i]
            p2 = nx_pos[j]
            if (i, j) in [(0, 12), (10, 5), (11, 6)]:
                edg4.append(DashedLine([p1, p2], spacing=0.3).lw(2))
            else:
                edg4.append(Line([p1, p2]).lw(2))

        for n in G.nodes():
            if n == 2:
                pts4.append(Sphere(r=r1[n]).pos(pos1[n]).alpha(0.4))
            elif n in [1, 3, 4]:
                pts4.append(Cylinder(r=r1[n], height=12.4, cap=True, pos=pos1[n], c=color[n]).alpha(0.2))

            if n in [0, 1, 2, 3, 4, 5, 6]:
                sph1.append(Sphere(r=0.5, pos=nx_pos[n], c='k'))
            else:
                sph1.append(Sphere(r=0, pos=nx_pos[n], c='k'))

            # sph1 = [Sphere(r=0.5, pos=nx_pos[n], c='k') for n in G.nodes]
        # axes = Axes(edg, xyGrid=True, zxGrid=True, yzGrid=True)
        plt.show(pts4, edg4, sph1, axes=False, interactive=True, bg='w', at=2)

        # --------------------------------------------------------------------------------------------------------------
        # subfigure 3
        # --------------------------------------------------------------------------------------------------------------
        pts1, edg1 = [], []
        for i, j in G.edges():
            p1 = nx_pos[i]
            p2 = nx_pos[j]
            if (i, j) in [(0, 12), (10, 5), (11, 6)]:
                edg1.append(DashedLine([p1, p2], spacing=0.3).lw(2))
            else:
                edg1.append(Line([p1, p2]).lw(2))

        for n in G.nodes():
            if n == 2:
                pts1.append(Sphere(r=r1[n]).pos(pos1[n]).alpha(0.4))
            elif n in [1, 3, 4]:
                cylb = Cylinder(r=r1[n], height=12.4, cap=True, pos=pos1[n], c=color[n]).alpha(0.2)
                cylt = Cylinder(r=r2[n], height=12.4, cap=True, pos=pos1[n], c=color[n]).alpha(0.2)
                cyl = merge(cylb, cylt)
                pts1.append(cyl)

            if n in [0, 1, 2, 3, 4, 5, 6]:
                sph1.append(Sphere(r=0.5, pos=nx_pos[n], c='k'))
            else:
                sph1.append(Sphere(r=0, pos=nx_pos[n], c='k'))

            # sph1 = [Sphere(r=0.5, pos=nx_pos[n], c='k') for n in G.nodes]
        # axes = Axes(edg, xyGrid=True, zxGrid=True, yzGrid=True)
        plt.show(pts1, edg1, sph1, axes=False, interactive=True, bg='w', at=3).screenshot("method")


    def Fig11(self):
        """ create interpolated volume for glc, lac concentration data
        Example:
        -------
        G = nx.gnm_random_graph(n=30, m=55, seed=1)
        nxpos = nx.spring_layout(G, dim=3)
        nxpts = [nxpos[pt] for pt in sorted(nxpos)]

        # node values
        vals = range(10, 10 + len(nxpts))
        print(vals)
        nx_lines, vals_segments = [], []
        for i, j in G.edges():
            nx_lines.append([nxpos[i], nxpos[j]])
            vals_segments += [vals[i], vals[j]]

            nx_pts = Points(nxpts, r=12)
            nx_edg = Lines(nx_lines).lw(2)

            nx_pts.cmap('YlGn', vals).addScalarBar()
            nx_edg.cmap('YlGn', vals_segments)
            labs = nx_pts.labels(vals, scale=0.03, precision=2).c('w')

        vol = interpolateToVolume(nx_pts, kernel='shepard', N=2, dims=(50, 50, 50))
        slices = []
        for i in range(0, 51, 10):
            print(i)
            sl = vol.xSlice(i).lighting('off').alpha(0.2).cmap('YlGn')
            slices.append(sl)
            slices.append(sl.isolines(vmin=10, vmax=40))
        show(nx_pts, nx_edg, labs, slices, bg='w', axes=9)
        # show(labs, slices, nx_pts.scalarbar, bg='w', axes=9, )
        """
        input = get_graph(sheet="test11") #sheet="test11_default"
        pos = dict(zip(input['nodes'], input['xyz']))
        lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
        pts = Points(input['xyz'], r=10).lighting('off').rotateX(-90)
        edg = Lines(lines).lw(5).rotateX(-90)

        mets = {"glc_ext": 0, "glc": 1, "lac": 3, "lac_ext": 2}
        color = {"glc_ext": 'Reds', "glc": 'Reds', "lac_ext": 'Blues', "lac": 'Blues'}

        data = {}
        time, scan = get_dynamic_data(fs=RESULTS_SIMGRAPH_FILE)
        bounds = graph_utils_py.get_bounds_pscan(scan=scan, common=True)
        print(bounds)
        x0, x1, y0, y1, z0, z1 = pts.bounds()
        for species, species_s in scan['simgraph'].items():
            data[species] = dict(zip(time, species_s))

        t = 30.0  #30.0
        species = "glc"
        txt = Text(
            f'[{species} ] (mM)  t = {t}s',
            font='VictorMono',
            justify='bottom-center',
            s=5,  # 40
            c='k'
        )
        txt.pos((x0 + x1) / 2, y0, z1)

        # --------------------------------------------------------------------------------------------------------------
        # cell species are plotted at points
        # --------------------------------------------------------------------------------------------------------------
        def cell_species(species, show=False):
            pts1 = pts.clone().cmap(color[species],
                                    data[species][t],
                                    vmin=bounds[species]['simgraph']['lb'],
                                    vmax=bounds[species]['simgraph']['ub'])
            pts1.addScalarBar3D(
                titleFont='VictorMono',
                labelFont='VictorMono',
                c='k'
            )
            pts1.scalarbar.addPos(0.5, 0, 0).useBounds()
            for a in pts1.scalarbar.unpack():
                if a and a.name == 'Text': a.scale(2.5)

            edg1 = edg.clone()

            # --------------------------------------------------------------------------------------------------------------
            #  method1
            if False:
                vol = interpolateToVolume(pts1, kernel='shepard', N=2, dims=(50, 50, 50))
                vol.c(color[species]).mode(0).alphaUnit(1).printInfo()
                if show: show(vol, pts1, edg1, bg='w', axes=9)
                return vol

            # method2
            if False:
                vol = interpolateToVolume(pts1, radius=5, kernel='linear', N=2)
                if species == "glc":
                    iso = vol.isosurface([2, 3, 4, 5]).alpha(0.1).cmap(color[species])  # [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5.0]
                elif species == "lac":
                    iso = vol.isosurface([1, 2, 3, 4]).alpha(0.2).cmap(color[species])
                if show: show(iso, pts1, edg1, bg='w', axes=True)
                return iso

            # method3
            if True:
                vol = interpolateToVolume(pts1, radius=5, kernel='linear', N=2, dims=(50, 50, 50))
                slices = []
                for i in range(0, 51, 5):
                    sl = vol.xSlice(i).lighting('off').alpha(0.5).cmap(color[species])
                    slices.append(sl)
                    # slices.append(sl.isolines(vmin=0, vmax=5))
                if show: show(pts1, edg1, slices, bg='w', axes=9)
                return slices

        # --------------------------------------------------------------------------------------------------------------
        # ext species are plotted at points and interpolated to segments
        # --------------------------------------------------------------------------------------------------------------

        def ext_species(species):
            pts1 = pts.clone().cmap(color[species],
                                    data[species][t],
                                    vmin=bounds[species]['simgraph']['lb'],
                                    vmax=bounds[species]['simgraph']['ub'])
            pts1.pointSize(0)
            pts1.addScalarBar3D(
                titleFont='VictorMono',
                labelFont='VictorMono',
                c='k'
            )
            pts1.scalarbar.addPos(0.1, 0, 0).useBounds()
            for a in pts1.scalarbar.unpack():
                if a and a.name == 'Text': a.scale(2.5)

            edg1 = edg.clone().cmap(color[species],
                                    get_vals_segment(point_data=data[species][t], sheet="test11"),
                                    vmin=bounds[species]['simgraph']['lb'],
                                    vmax=bounds[species]['simgraph']['ub'])
            return pts1, edg1
        # show(pts1, edg1, axes=True, bg2='white', at=0, interactive=True)

        vol_glc = cell_species(species="glc")
        pts_eglc, edg_eglc = ext_species(species="glc_ext")
        vol_lac = cell_species(species="lac")
        pts_elac, edg_elac = ext_species(species="lac_ext")

        plt12 = Plotter(N=2, interactive=True)
        plt12.show(pts_eglc, edg_eglc, vol_glc, at=0) #vol_glc
        plt12.show(pts_elac, edg_elac, vol_lac, at=1).screenshot("vol2") #vol_lac

    def Fig12(self):
        """
        :param self:
        :return:
        reference:
        ---------
        https://github.com/marcomusy/vedo/blob/master/examples/volumetric/probeLine2.py
        """
        input = get_graph(sheet="test11")
        pos = dict(zip(input['nodes'], input['xyz']))
        lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
        pts = Points(input['xyz'], r=10).lighting('off')
        edg = Lines(lines).lw(5)

        color = {"glc_ext": 'Reds', "glc": 'Reds', "lac_ext": 'Blues', "lac": 'Blues'}
        time, scan = get_dynamic_data(fs=RESULTS_SIMGRAPH_FILE)
        bounds = graph_utils_py.get_bounds_pscan(scan=scan, common=True)

        data = {}
        for species, species_s in scan['simgraph'].items():
            data[species] = dict(zip(time, species_s))
        ts = [0.0, 1.2, 4.2, 15, 30, 60, 120, 300]

        y = {}
        for t in ts:
            species = "glc"
            pts1 = pts.clone().cmap(color[species],
                                    data[species][t],
                                    vmin=bounds[species]['simgraph']['lb'],
                                    vmax=bounds[species]['simgraph']['ub'])
            pts1.addScalarBar3D(
                titleFont='VictorMono',
                labelFont='VictorMono',
                c='k'
            )
            pts1.scalarbar.addPos(0.1, 0, 0).useBounds()
            for a in pts1.scalarbar.unpack():
                if a and a.name == 'Text': a.scale(2.5)

            vol = interpolateToVolume(pts1, radius=5, kernel='linear', N=2, dims=(50, 50, 50))
            p1, p2 = pos[35], pos[38]
            pl = probeLine(vol, p1, p2, res=10).lineWidth(4).color('r')

            x = pl.points()
            y.update({t: pl.getPointArray()})

        ax = plt.subplot(111)
        for t in ts:
            ax.plot(
                x[:, 0], y[t],
                # color=color[n],
                markersize=5,
            )
        ax.set_xlabel('Time (s)', fontsize=10)
        ax.set_ylabel('Concentration (mM)', fontsize=10)
        plt.show()

        show(pts, edg, pl, bg='w')

    def Fig13(self):
        """plots steady-state time  as node values and interpolates"""

        def plot_risetime(key):
            data = scan[key]
            node_attr, edge_attr = data['node_attr'], data['edge_attr']

            vplt = Plotter(interactive=False)
            # vmin, vmax = get_bounds_static(scan=scan, param='tss_glc_ext', attr='node_attr')

            pts1 = pts.clone().pointSize(10).pointColors(node_attr['tss_glc_ext'], cmap='Spectral')
            pts1.addScalarBar3D(
                title='(a) Rise time (s)',
                titleFont='VictorMono',
                labelFont='VictorMono',
                c='k',
            )
            pts1.scalarbar.addPos(1, 0, 0).useBounds()

            label = [str(n) for n in input['nodes']]
            labs = pts.labels(label, scale=5, font='VictorMono', c='k').addPos(0.05, 0, 0.5)

            edg1 = edg.clone().pointColors(edge_attr['tss_glc_ext'], cmap='Spectral')
            for a in pts1.scalarbar.unpack():
                if a and a.name == 'Text': a.scale(1.5)

            if test == "test11_default":
                axes = Axes(edg1, xtitle='y', ytitle='z', ztitle='x', xyGrid=True, zxGrid=True, yzGrid=True)
            else:
                axes = True

            vplt_return = vplt.show(pts1, edg1, labs, axes)
            return vplt_return

        input = get_graph(sheet="test24")
        # create vtk polydata
        pos = dict(zip(input['nodes'], input['xyz']))
        lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
        pts = Points(input['xyz'], r=10).lighting('off').rotateX(0)  #-90 for test11
        edg = Lines(lines).lw(5).rotateX(0)

        scan = get_static_data(fs=RESULTS_SIMGRAPH_FILE, fc=RESULTS_COMSOL_FILE, twitch=True)
        vplt_sim = plot_risetime(key="simgraph")
        vplt_com = plot_risetime(key="comsol")

        f = os.path.join(RESULTS_PLOTS_DIR, f'{test}_risetime.png')
        plt12 = Plotter(N=2, interactive=True, sharecam=True)
        plt12.show(vplt_sim.actors, at=0, zoom=1.2)
        plt12.show(vplt_com.actors, at=1, zoom=1.2).screenshot(f)

    def Fig14(self):
        """plots rise time  as node values and interpolates"""

        def risetime(sheet):
            input = get_graph(sheet=sheet)
            pos = dict(zip(input['nodes'], input['xyz']))
            lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
            pts = Points(input['xyz'], r=10).lighting('off').rotateX(0)
            edg = Lines(lines).lw(5).rotateX(0)

            data = scan[sheet]
            node_attr, edge_attr = data['node_attr'], data['edge_attr']

            pts1 = pts.clone().pointSize(0).cmap(lut, node_attr['rt_glc_ext'])
            if sheet == "test9":
                pts1.addScalarBar3D(
                    title='Rise time (s)',
                    titleFont='VictorMono',
                    labelFont='VictorMono',
                    c='k',
                )
                pts1.scalarbar.addPos(40, 0, 0).useBounds()
                for a in pts1.scalarbar.unpack():
                    if a and a.name == 'Text': a.scale(1.5)

            label = [str(n) for n in input['nodes']]
            labs = pts.labels(label, scale=5, font='VictorMono', c='k').addPos(0.05, 0, 0.5)

            edg1 = edg.clone().cmap(lut, edge_attr['rt_glc_ext'])

            f = os.path.join(RESULTS_PLOTS_DIR, f'{sheet}_risetime.png')
            vplt = Plotter(interactive=True)
            vplt_return = vplt.show(pts1, edg1, axes=True, zoom=zoom[sheet]).screenshot(f)
            return vplt_return

        scan = {}
        tag = {"test10": "design2", "test9": "design1"}
        zoom = {"test10": 1, "test9": 1}
        for s in tag.keys():
            fs = os.path.join(RESULTS_SIMGRAPH_DIR, f"{s}_default_bc1_v1_c0.mat")
            r = get_static_data(fs=fs, twitch=True, sheet=s)
            scan[s] = r["simgraph"]

        vmin, vmax = get_bounds_static(scan=scan, param='rt_glc_ext', attr='node_attr')  # vmin=0, vmax=2034.1

        # manually build a lookup table for mapping colors
        lut = buildLUT([
            (vmin, 'db'),
            (15, 'b'),
            (35, 'lb'),
            (54, 'w'),
            (vmax, 'r'),
        ],
            vmin=vmin, vmax=vmax,
            interpolate=True
        )

        plt1 = risetime(sheet="test10")
        plt2 = risetime(sheet="test9")
        fpath = os.path.join(RESULTS_PLOTS_DIR, 'risetime.png')
        plt12 = Plotter(N=2, interactive=True, sharecam=False)
        plt12.show(plt1.actors, at=0, zoom=1, axes=True)
        plt12.show(plt2.actors, at=1, zoom=1, axes=True).screenshot(fpath)

    def Fig15(self):
        """ creates time-series plots of blood and cell species for islet vasculature"""
        settings.showRendererFrame = False

        input = get_graph(sheet="test11")
        nodes = input['nodes']
        pos = dict(zip(nodes, input['xyz']))
        lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
        pts = Points(input['xyz'], r=10).lighting('off').rotateX(-90)
        edg = Lines(lines).lw(5).rotateX(-90)
        rxn_nodes = [1, 11, 15, 17, 18, 24, 28, 30, 33, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
                     52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
                     74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                     96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112]
        sph_r = OrderedDict.fromkeys(nodes, 0)
        sph_r.update(OrderedDict.fromkeys(rxn_nodes, 3))

        color = {"glc_ext": 'Reds', "glc": 'Reds', "lac_ext": 'Blues', "lac": 'Blues'}
        label = {"glc_ext": 'eglucose', "glc": 'glucose',
                 "lac_ext": 'elactate', "lac": 'lactate'}
        data = {}
        time, scan = get_dynamic_data(fs=RESULTS_SIMGRAPH_FILE)
        ts = [0.6, 6.0, 30.0]
        plt = Plotter(shape=(len(ts), nspecies), interactive=0, sharecam=1, bg2='white')  #

        # ------------------------------------------------------------------------------------------------------------------
        x0, x1, y0, y1, z0, z1 = pts.bounds()
        bounds = graph_utils_py.get_bounds_pscan(scan=scan, common=False)

        at = 0
        for idx, key in enumerate(scan):
            for i in range(nspecies):
                s = SPECIES[i]
                c = scan[key][s]
                data[s] = dict(zip(time, c))

            for t in ts:
                for species in SPECIES:
                    txt = Text(
                        f't: {t}s',
                        font='VictorMono',
                        justify='bottom-left',
                        s=5,
                        c='k'
                    )
                    txt.pos((x0 + x1) / 2, y0, z1)
                    vmin = bounds[species]['simgraph']['lb']
                    vmax = bounds[species]['simgraph']['ub']

                    if species in ["glc", "lac"]:
                        temp = []
                        for value in data[species][t]:
                            value_c = colorMap(
                                value,
                                name=color[species],
                                vmin=vmin,
                                vmax=vmax
                            )
                            temp.append(value_c)
                        sph_c = OrderedDict(zip(nodes, temp))

                        # pts
                        pts1 = pts.clone().cmap(
                            color[species],
                            data[species][t],
                            vmin=vmin,
                            vmax=vmax
                        )

                        sph1 = [Sphere(r=sph_r[n]).color(sph_c[n]).pos(pos[n]).alpha(1).rotateX(-90) for n in nodes]
                        if ts.index(t) <= len(ts)-1:
                            pts1.addScalarBar3D(
                                title=f'[{label[species]}] (mM)',
                                titleFont='VictorMono',
                                labelFont='VictorMono',
                                c='k'
                            )
                            pts1.scalarbar.addPos(15, 0, 0).useBounds()
                            for a in pts1.scalarbar.unpack():
                                if a and a.name == 'Text': a.scale(2.5)

                        # edg
                        edg1 = edg.clone()

                    else:
                        pts1 = pts.clone().pointSize(2).cmap(
                            color[species],
                            data[species][t],
                            vmin=vmin,
                            vmax=vmax
                        )
                        edg1 = edg.clone().cmap(
                            color[species],
                            get_vals_segment(point_data=data[species][t], sheet="test11"),
                            vmin=vmin,
                            vmax=vmax
                        )

                        if ts.index(t) <= len(ts)-1:
                            edg1.addScalarBar3D(
                                title=f'[{label[species]}] (mM)',
                                titleFont='VictorMono',
                                labelFont='VictorMono',
                                c='k'
                            )
                            edg1.scalarbar.addPos(15, 0, 0).useBounds()
                            for a in edg1.scalarbar.unpack():
                                if a and a.name == 'Text':
                                    a.scale(2)

                    axes = Axes(edg1, xtitle='y', ytitle='z', ztitle='x', xyGrid=True, zxGrid=False, yzGrid=False)
                    if species in ["glc", "lac"]:
                        plt.show(
                            # pts1,
                            sph1,
                            edg1,
                            pts1.scalarbar,
                            txt,
                            axes,
                            at=at,
                            zoom=1,
                        )
                    else:
                        plt.show(
                            pts1,
                            edg1,
                            txt,
                            axes,
                            at=at,
                            zoom=1,
                        )

                    at = at + 1
            interactive()
            plt.screenshot(os.path.join(RESULTS_PLOTS_DIR, f'{test}_conc_bvcell.png'))

    def Fig16(self):
        """
        plot glucose and lactate net exchange flux
        Reference:
        https://stackoverflow.com/questions/14762181/adding-a-y-axis-label-to-secondary-y-axis-in-matplotlib
        https://stackoverflow.com/questions/58490203/matplotlib-plot-a-line-with-open-markers-where-the-line-is-not-visible-within
        ---------
        Example:
        x = np.arange(0, 10, 0.1)
        y1 = 0.05 * x ** 2
        y2 = -1 * y1
        """

        data = io_utils_py._get_file(file=RESULTS_GSCAN_FILE, spreadsheet=False)

        dose, glcim, lacex = [], [], []
        for s in data['variable']:
            dose.append(s.dose)
            glcim.append(s.glcim_net)
            lacex.append(abs(s.lacex_net))

        fig, ax = plt.subplots()
        # fig, ax1 = plt.subplots()
        # ax2 = ax1.twinx()
        ax.plot(dose, glcim, 'r', marker='o', markerfacecolor='none', markersize=4)
        ax.plot(dose, lacex, 'b', marker='s', markerfacecolor='none', markersize=4)

        ax.set_xlabel('Glucose dose (mM)', fontsize=12)
        ax.set_ylabel('Net exchange flux (pmole/min)', fontsize=12)
        # ax2.set_ylabel('Lactate release (pmole/min)', fontsize=12)

        fig.show()
        f = os.path.join(RESULTS_PLOTS_DIR, f"{test}_flux.svg")
        ax.figure.savefig(
            f,
            transparent=True,
            dpi=600,
            bbox_inches="tight"
        )

    def Fig17(self):
        """
        compare conc simgraph vs. comsol for glc_ext evolution in islet vasculature
        :param pts:
        :param edg:
        :return:
        ref: https://stackoverflow.com/questions/44947505/how-to-make-a-movie-out-of-images-in-python
        Note to self: If the colormap appears black/blank for cell species, check input file
        """

        def plot_timeseries(scan):

            images = []
            data = {}
            for key, species_s in scan.items():
                data_s = dict(zip(time, species_s["glc_ext"]))
                data[key] = data_s

            # ts = [0.0, 1.2, 4.2]
            bounds = graph_utils_py.get_bounds_pscan(scan=scan, common=True)
            x0, x1, y0, y1, z0, z1 = pts.bounds()
            plt = Plotter(N=2, interactive=False)

            for t in time:
                for key in scan.keys():
                    txt = Text3D(
                        f'[glc_ext ] (mM)  t = {round(t,2)}s',
                        font='VictorMono',
                        justify='bottom-center',
                        s=5,
                        c='k'
                    )
                    txt.pos((x0 + x1) / 2, y0, z1)

                    # pts
                    pts1 = pts.clone()
                    pts1.pointSize(1)

                    # edg
                    edg1 = edg.clone().cmap(
                        "Reds",
                        get_vals_segment(point_data=data[key][t], sheet="test11"),
                        vmin=bounds["glc_ext"][key]['lb'],
                        vmax=bounds["glc_ext"][key]['ub']
                    )
                    edg1.addScalarBar3D(
                        title=key,
                        titleFont='VictorMono',
                        labelFont='VictorMono',
                        c='k',
                        titleXOffset=3.5
                    )
                    edg1.scalarbar.addPos(0.1, 0, 0).useBounds()
                    for a in edg1.scalarbar.unpack():
                        if a and a.name == 'Text': a.scale(2.5)

                    if key == "simgraph":
                        plt.show(
                            pts1,
                            edg1,
                            txt,
                            axes=True,
                            bg2='white',
                            at=0,
                            # zoom=1.5,
                        )
                    elif key == "comsol":
                        plt.show(
                            pts1,
                            edg1,
                            txt,
                            axes=True,
                            bg2='white',
                            at=1,
                            # zoom=1.5,
                        )

                    fpath = os.path.join(RESULTS_PLOTS_DIR, f'{test}_conc_glc_ext_{round(t,2)}s.png')
                images.append(fpath)
                plt.screenshot(fpath)

            return images

        input = get_graph(sheet="test11")
        pos = dict(zip(input['nodes'], input['xyz']))
        lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
        pts = Points(input['xyz'], r=10).lighting('off').rotateX(-90)
        edg = Lines(lines).lw(5).rotateX(-90)

        time, scan = get_dynamic_data(fs=RESULTS_SIMGRAPH_FILE, fc=RESULTS_COMSOL_FILE)

        # create video
        images = plot_timeseries(scan=scan)
        print(images)
        video_name = os.path.join(RESULTS_PLOTS_DIR, "islet_compare.mp4")
        frame = cv2.imread(os.path.join(RESULTS_PLOTS_DIR, images[0]))
        height, width, layers = frame.shape
        fourcc = 0
        fps = 1
        video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))
        for image in images:
            fpath = os.path.join(RESULTS_PLOTS_DIR, image)
            video.write(cv2.imread(fpath))
            os.remove(fpath)  # remove the image file
        cv2.destroyAllWindows()
        video.release()

    def Fig18(self):
        """
        compare conc simgraph vs. comsol for glc_ext evolution in mesentery vasculature
        :param pts:
        :param edg:
        :return:
        ref: https://stackoverflow.com/questions/44947505/how-to-make-a-movie-out-of-images-in-python
        Note to self: If the colormap appears black/blank for cell species, check input file
        """

        def plot_timeseries(scan, sheet):

            images = []
            data = {}
            for key, species_s in scan.items():
                data_s = dict(zip(time, species_s["glc_ext"]))
                data[key] = data_s

            # ts = [0.0, 1.2, 4.2]
            bounds = graph_utils_py.get_bounds_pscan(scan=scan, common=True)
            x0, x1, y0, y1, z0, z1 = pts.bounds()
            plt = Plotter(N=2, interactive=False)


            for t in time:
                for key in scan.keys():
                    txt = Text3D(
                        f'[glc_ext ] (mM)  t = {round(t,2)}s',
                        font='VictorMono',
                        justify='bottom-center',
                        s=150,
                        c='k'
                    )
                    txt.pos((x0 + x1) / 2, y0, 5000)

                    # pts
                    pts1 = pts.clone()
                    pts1.pointSize(1)

                    # edg
                    edg1 = edg.clone().cmap(
                        lut,
                        get_vals_segment(point_data=data[key][t], sheet=sheet),
                        vmin=bounds["glc_ext"][key]['lb'],
                        vmax=bounds["glc_ext"][key]['ub']
                    )
                    edg1.addScalarBar3D(
                        title=key,
                        titleFont='VictorMono',
                        labelFont='VictorMono',
                        c='k',
                        titleXOffset=2.5
                    )
                    edg1.scalarbar.addPos(0.1, 0, 0).useBounds()
                    for a in edg1.scalarbar.unpack():
                        if a and a.name == 'Text': a.scale(1)

                    if key == "simgraph":
                        plt.show(
                            pts1,
                            edg1,
                            txt,
                            axes=True,
                            bg2='white',
                            at=0,
                            zoom=1.3,
                        )
                    elif key == "comsol":
                        plt.show(
                            pts1,
                            edg1,
                            txt,
                            axes=True,
                            bg2='white',
                            at=1,
                            zoom=1.3,
                        )

                    fpath = os.path.join(RESULTS_PLOTS_DIR, f'{test}_conc_glc_ext_{round(t,2)}s.png')
                images.append(fpath)
                plt.screenshot(fpath)

            return images

        input = get_graph(sheet="test24")
        pos = dict(zip(input['nodes'], input['xyz']))
        lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
        pts = Points(input['xyz'], r=10).lighting('off')
        edg = Lines(lines).lw(5)

        # ------------------------------------------------------------------------------------------------------------------
        # manually build a lookup table for mapping colors
        # ------------------------------------------------------------------------------------------------------------------
        lut = buildLUT([
            (0, 'w'),
            (5, 'red'),
        ],
            vmin=0,
            vmax=5,
            interpolate=True,
            aboveColor='darkred',
            belowColor='w'
        )
        # ------------------------------------------------------------------------------------------------------------------

        time, scan = get_dynamic_data(fs=RESULTS_SIMGRAPH_FILE, fc=RESULTS_COMSOL_FILE)

        # create video
        images = plot_timeseries(scan=scan, sheet="test24")
        video_name = os.path.join(RESULTS_PLOTS_DIR, "mesentery_compare.mp4")
        frame = cv2.imread(os.path.join(RESULTS_PLOTS_DIR, images[0]))
        height, width, layers = frame.shape
        fourcc = 0
        fps = 1
        video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))
        for image in images:
            fpath = os.path.join(RESULTS_PLOTS_DIR, image)
            video.write(cv2.imread(fpath))
            os.remove(fpath)  # remove the image file
        cv2.destroyAllWindows()
        video.release()


    def Fig19(self):
        """
        compare conc simgraph vs. comsol for glc_ext evolution in tumor vasculature
        :param pts:
        :param edg:
        :return:
        ref: https://stackoverflow.com/questions/44947505/how-to-make-a-movie-out-of-images-in-python
        Note to self: If the colormap appears black/blank for cell species, check input file
        """

        def plot_timeseries(scan, sheet, t):

            input = get_graph(sheet=sheet)
            pos = dict(zip(input['nodes'], input['xyz']))
            lines = [(pos[t], pos[h]) for t, h in input['ed_ls']]
            pts = Points(input['xyz'], r=10).lighting('off').rotateX(0)
            edg = Lines(lines).lw(5).rotateX(0)

            data = dict(zip(time, scan["simgraph"]["glc_ext"]))

            bounds = graph_utils_py.get_bounds_pscan(scan=scan, common=True)
            x0, x1, y0, y1, z0, z1 = pts.bounds()
            plt = Plotter(N=2, interactive=False)
            label = {"test9": "Tumor design1", "test10": "Tumor design2"}

            for key in scan.keys():
                txt = Text3D(
                    f'[glc_ext ] (mM)  t = {round(t,2)}s',
                    font='VictorMono',
                    justify='bottom-center',
                    s=25,
                    c='k'
                )
                txt.pos((x0 + x1) / 2, y0, z1+500)

                # pts
                pts1 = pts.clone()
                pts1 = pts.clone().cmap(
                    "Reds",
                    data[t],
                    vmin=bounds["glc_ext"][key]['lb'],
                    vmax=bounds["glc_ext"][key]['ub']
                )
                pts1.pointSize(1)
                # edg
                edg1 = edg.clone().cmap(
                    "Reds",
                    get_vals_segment(point_data=data[t], sheet=sheet),
                    vmin=bounds["glc_ext"][key]['lb'],
                    vmax=bounds["glc_ext"][key]['ub']
                )
                edg1.addScalarBar3D(
                    title=label[sheet],
                    titleFont='VictorMono',
                    labelFont='VictorMono',
                    c='k',
                    titleXOffset=2.5
                )

                edg1.scalarbar.addPos(0.1, 0, 0).useBounds()
                for a in edg1.scalarbar.unpack():
                    if a and a.name == 'Text': a.scale(1.5)

                vplt_return = plt.add(
                    [pts1,
                    edg1,
                    txt],
                    # axes=True,
                    # bg2='white',
                    at=0,
                    # zoom=1.5,
                )

            return vplt_return

        scan = {}
        time, scan["design1"] = get_dynamic_data(fs=os.path.join(RESULTS_SIMGRAPH_DIR, "test9_default_bc1_v1_c0.mat"))
        time, scan["design2"] = get_dynamic_data(fs=os.path.join(RESULTS_SIMGRAPH_DIR, "test10_default_bc1_v1_c0.mat"))
        images = []
        for t in time:
            plt1 = plot_timeseries(scan=scan["design1"], sheet="test9", t=t)
            plt2 = plot_timeseries(scan=scan["design2"], sheet="test10",  t=t)
            fpath = os.path.join(RESULTS_PLOTS_DIR, f'tumor_conc_glc_ext_{round(t,2)}s.png')
            plt12 = Plotter(N=2, interactive=False, sharecam=True)
            plt12.show(plt1.actors, at=0, zoom=1.2, axes=True, bg2='gray')
            plt12.show(plt2.actors, at=1, zoom=1.2, axes=True, bg2='gray').screenshot(fpath)
            images.append(fpath)

        # create video
        video_name = os.path.join(RESULTS_PLOTS_DIR, "tumor_compare.mp4")
        frame = cv2.imread(os.path.join(RESULTS_PLOTS_DIR, images[0]))
        height, width, layers = frame.shape
        fourcc = 0
        fps = 1
        video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))
        for image in images:
            fpath = os.path.join(RESULTS_PLOTS_DIR, image)
            video.write(cv2.imread(fpath))
            os.remove(fpath)  # remove the image file
        cv2.destroyAllWindows()
        video.release()


if __name__ == '__main__':

    # ------------------------------------------------------------------------------------------------------------------
    # Generate geometry
    # ------------------------------------------------------------------------------------------------------------------
    G = create_graph(directed=True, attributes=True)
    # G = get_eucledian_lengths_wrt_origin(origin=34, edge_ed=False)
    # G, H = mesh_network(G=G)
    # draw_graph3d_vedo(Gs=[G], label_scale=5, highlight=False, default_label=True)
    # exit()
    manuscript_figures = ManuscriptFigures()
    manuscript_figures.figures()
