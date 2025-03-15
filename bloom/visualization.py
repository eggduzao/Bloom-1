"""
visualization
-------------

# Example Usage
vp = VizPlots()
vp.pie_chart(values=[10, 20, 30], labels=['A', 'B', 'C'])
vp.star_plot({'Metric1': 3, 'Metric2': 5, 'Metric3': 4})
vp.bar_plot([5, 6, 7], x_labels=['X', 'Y', 'Z'])
vp.box_plot([[1, 2, 3], [3, 4, 5]])
vp.distribution_plot([[1, 2, 3, 2, 2], [3, 4, 5, 3, 3]])
vp.polar_plot({'Metric1': 1, 'Metric2': 4, 'Metric3': 3})
"""

###############################################################################
# Imports
###############################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.special import expit, softmax
import seaborn as sns
from collections import OrderedDict

###############################################################################
# Constants
###############################################################################

# Constants
SEED = 1987
np.random.seed(SEED)

###############################################################################
# Classes
###############################################################################

class MetricTables:
    """
    Metric                                          A. Range    B. Best?    C. Per Image / Tile / Cell?                         D. Best Plot Type
    1. Dice Similarity Coefficient (DSC)            [0, 1]      Increase    Per structure (cell/nucleus/tissue region)          Violin (distribution)
    2. Intersection over Union (IoU, Jaccard Index) [0, 1]      Increase    Per structure (cell/nucleus/tissue region)          Violin (distribution)
    3. Pixel Accuracy (PxAcc)                       [0, 1]      Increase    Per image/tile (DAB, SegPath, disease prediction)   Violin (distribution) or Bar (mean)
    4. Hausdorff Distance (HD)                      [0, ∞]      Decrease    Per structure (cell/nucleus/tissue region)          Violin (distribution)
    5. Adjusted Rand Index (ARI)                    [-1, 1]     Increase    Per structure (cell/nucleus segmentation)           Violin (distribution)
    ---
    Metric                                      A. Range        B. Best?    C. Per Image / Tile / Run?          D. Best Plot Type
    6. Inference Time (Seconds per Image)       [0, ∞]          Decrease    Per image or per batch (run)        Bar (mean + error bars)
    7. Memory Usage (RAM & VRAM Consumption)    [0, ∞] (GB)     Decrease    Per run (entire process)            Line (vs. batch size) or Bar
    8. Model Size (MB)                          [0, ∞] (MB)     Decrease    Per model (global, not per image)   Bar (sorted by size)
    9. Energy Consumption (Watts per Image)     [0, ∞] (W)      Decrease    Per image or per batch (run)        Bar (mean + error bars)
    -----------
    SegPath:            Total: 292,398                  Train: 266,724                      Test: 25,674   
    SNOW:               Total: 20,000 tiles (512×512)   Total: 1,448,522 annotated nuclei
    NeurIPS Challenge:  Total: 1,659 Images
    TissueNet:          Total: 7,022 Images             Train: 2,580 512x512        Validation: 3,118 256x256   Test: 1,324 256x256
    DynamicNuclearNet Segmentation: Total: 7,084        Annotations: 913,119
    BCData:             Total: 1,338 Images
    PanNuke:            Total: 7,901 Images             Total: 189,744 annotated nuclei
    DAB Neila:          Total = 4,977                   Total PE: 2,004 / Total IMIP: 2,973 
    """

    def __init__(self):

        # Methods
        self.wyw_methods = ["GhostNet", "Cellpose2", "nnU_Net", "KIT_GE", "Omnipose", "Mesmer", "QuPath", "CelloType"]
        self.bio_methods ["SNOW", "CellProfler4", "HistomicsTK", "ImageJ", "ilastik", "EpidermaQuant", "DAB_quant"]
        self.nipschal_methods = ["osilab", "sribdmed", "cells", "saltfish", "redcat_autox", "train4ever", "overoverfitting", 
                                 "vipa", "naf", "bupt_mcprl", "cphitsz", "wonderworker", "cvmli", "m1n9x", "fzu312", "sgroup",
                                 "smf", "sanmed_ai", "hilab", "guanlab", "daschlab", "mbzuai_cellseg", "quiil", "plf", 
                                 "siatcct", "nonozz", "littlefatfish", "boe_aiot_cto"]
        self.comp_methods = ["SSTBM", "UNet", "Att_UNet", "SAP_UNet", "LUNet", "PsdLabeling", "NucCell_GAN", "DeepLabV3P", 
                             "HoVer_Net", "Swin", "Swin_V2", "ViT_MoE", "Random_Forest", "SVM", "Watershed", "Otsu"]
        self.all_methods = self.wyw_methods + self.bio_methods + self.nipschal_methods + self.comp_methods
        self.selected_methods = []

        # Datasets
        self.datasets = {"SegPath": 25_674, "SNOW": 20_000, "NeurIPS": 1_659, "TissueNet": 7_022, "DynamicNuclearNet": 7_084,
                         "BCData": 1_338, "PanNuke": 7_901, "DAB": 4_977}

        # Metrics
        self.dice = OrderedDict()
        self.iou = OrderedDict()
        self.pxAcc = OrderedDict()
        self.hd = OrderedDict()
        self.ari = OrderedDict()
        self.time = OrderedDict()
        self.mem = OrderedDict()
        self.modsize = OrderedDict()
        self.energy = OrderedDict()

    def create_segpath_table(self):
        """
        59 methods x 9 metrics x 25,674 images
        """

    def create_snow_table(self):

    def create_neurips_table(self):

    def create_tissuenet_table(self):

    def create_nuclearnet_table(self):

    def create_bcdata_table(self):

    def create_pannuke_table(self):

    def create_dab_table(self):


class VizPlots:
    """
    A comprehensive visualization class for generating various plots based on quality metrics.

    Methods
    -------
    pie_chart(data, labels, is_percent=False, plot_as_percent=False)
        Generates a pie chart with customizable formatting.
    star_plot(data)
        Generates a star plot with values from 1 to 5.
    bar_plot(data, x_labels, y_labels, x_title, y_title, inverted=False)
        Generates a bar plot with customizable aesthetics.
    box_plot(data)
        Generates a violin plot to represent distributions.
    distribution_plot(data)
        Generates overlapping distribution plots.
    polar_plot(data)
        Generates a polar plot.
    """

    def __init__(self):
        pass

    def pie_chart(self, values=None, labels=None, data=None, is_percent=False, plot_as_percent=False):
        """
        Generates a pie chart with optional data transformation to percentages.

        Parameters
        ----------
        values : list, optional
            List of numerical values.
        labels : list, optional
            List of labels for the pie chart.
        data : dict, optional
            Dictionary with labels as keys and numerical values as values.
        is_percent : bool, optional
            If True, data is already in percentage.
        plot_as_percent : bool, optional
            If True, display percentages on the pie chart.
        """
        if data:
            labels, values = zip(*data.items())

        if not is_percent:
            total = sum(values)
            values = [v / total * 100 for v in values]

        colors = sns.color_palette('pastel')
        explode = [0.05] * len(values)

        fig, ax = plt.subplots()
        wedges, texts, autotexts = ax.pie(
            values, labels=labels, autopct='%1.1f%%' if plot_as_percent else None,
            startangle=140, colors=colors, explode=explode, shadow=True
        )
        plt.title("Pie Chart")
        plt.show()

    def star_plot(self, data):
        """
        Generates a star plot for 1-to-5 ratings.

        Parameters
        ----------
        data : dict
            Dictionary with categories as keys and integer ratings (1-5) as values.
        """
        labels, values = zip(*data.items())

        fig, ax = plt.subplots(figsize=(6, 6))
        angles = np.linspace(0, 2 * np.pi, len(values), endpoint=False).tolist()
        values += (values[0],)
        angles += angles[:1]

        ax = plt.subplot(111, polar=True)
        ax.fill(angles, values, color='lightblue', alpha=0.6)
        ax.plot(angles, values, color='blue', linewidth=2)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels)
        plt.title("Star Plot")
        plt.show()

    def bar_plot(self, data, x_labels=None, y_labels=None, x_title="", y_title="", inverted=False):
        """
        Generates a bar plot with customizable options.

        Parameters
        ----------
        data : list
            List of numerical values or list of lists for series.
        x_labels : list, optional
            Labels for the x-axis.
        y_labels : list, optional
            Labels for the y-axis.
        x_title : str, optional
            Title for the x-axis.
        y_title : str, optional
            Title for the y-axis.
        inverted : bool, optional
            If True, inverts the axes.
        """
        if isinstance(data[0], (list, np.ndarray)):
            df = pd.DataFrame(data).T
            df.plot(kind="bar", stacked=True)
        else:
            plt.bar(x_labels, data)

        plt.xlabel(x_title)
        plt.ylabel(y_title)
        if inverted:
            plt.gca().invert_yaxis()
        plt.title("Bar Plot")
        plt.show()

    def box_plot(self, data, violin=True):
        """
        Generates box or violin plots for distributions.

        Parameters
        ----------
        data : list of lists
            List of distributions.
        violin : bool, optional
            If True, use violin plot instead of box plot.
        """
        df = pd.DataFrame(data).T
        if violin:
            sns.violinplot(data=df, palette='pastel')
        else:
            sns.boxplot(data=df, palette='pastel')
        plt.title("Box/Violin Plot")
        plt.show()

    def distribution_plot(self, data):
        """
        Generates distribution plots with overlapping histograms.

        Parameters
        ----------
        data : list of lists
            List of distributions.
        """
        df = pd.DataFrame(data).T
        for col in df.columns:
            sns.kdeplot(df[col], fill=True)
        plt.title("Distribution Plot")
        plt.show()

    def polar_plot(self, data):
        """
        Generates a polar plot for metrics.

        Parameters
        ----------
        data : dict
            Dictionary with metric names as keys and numerical values as values.
        """
        labels, values = zip(*data.items())
        angles = np.linspace(0, 2 * np.pi, len(values), endpoint=False).tolist()
        values += (values[0],)
        angles += angles[:1]

        ax = plt.subplot(111, polar=True)
        ax.fill(angles, values, color='lightgreen', alpha=0.6)
        ax.plot(angles, values, color='green', linewidth=2)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels)
        plt.title("Polar Plot")
        plt.show()

