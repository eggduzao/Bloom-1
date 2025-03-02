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
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

###############################################################################
# Constants
###############################################################################

# Constants
SEED = 1987
np.random.seed(SEED)

###############################################################################
# Classes
###############################################################################

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

