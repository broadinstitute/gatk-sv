import numpy
import scipy.integrate
import bisect
from tqdm import tqdm
import operator
import matplotlib
import seaborn
from matplotlib import pyplot
from matplotlib import colors
from matplotlib import patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import pandas
from typing import Optional, Union, Tuple, List, Sequence, Mapping, Text, Any
from types import MethodType, MappingProxyType


Numeric = Union[int, float, numpy.integer, numpy.floating]
Record = Union[numpy.record, Tuple]
Vector = Union[pandas.Series, numpy.ndarray]
ListLike = Union[List, pandas.Series, numpy.ndarray]


class Default:
    tick_labelrotation = -30
    num_1d_bins = 100
    num_2d_bins = 50
    min_1d_hist_width = 1.0e-6
    title_fontsize = 16
    label_fontsize = 16
    num_palette_colors = 12
    color_palette = "colorblind"
    # Default params overrides for matplotlib
    matplotlib_params = MappingProxyType({
        'savefig.dpi': 150,  # to adjust notebook inline plot size
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'font.size': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.figsize': [9.75, 6]
    })
    use_seaborn = True


def init_plotting_environment(
        matplotlib_params: Mapping[Text, Any] = Default.matplotlib_params,
        use_seaborn: bool = Default.use_seaborn,
        color_palette: str = Default.color_palette,
        num_palette_colors: int = Default.num_palette_colors
):
    # restore defaults if you've been mucking around
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    if use_seaborn:
        seaborn.set_palette(palette=color_palette, n_colors=num_palette_colors)
        seaborn.set(color_codes=True)
        pyplot.style.use('seaborn')
    matplotlib.rcParams.update(matplotlib_params)


init_plotting_environment()


def next_fig() -> int:
    """
    Get appropriate handle for next figure
    """
    fignums = pyplot.get_fignums()
    if fignums:
        fignums = sorted(fignums)
        n_figs = len(fignums)
        if fignums[-1] == n_figs:
            next_num = n_figs + 1
        else:
            next_num = -1
            for n, fignum in enumerate(fignums, start=1):
                if fignum > n:
                    next_num = n
                    break
            if next_num <= 0:
                raise RuntimeError("Unable to find next figure number")
    else:
        next_num = 1
    return next_num


def ordered_legend(ax: pyplot.Axes, labels_in_order: Sequence, **legend_kwargs):
    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: labels_in_order.index(t[0])))
    ax.legend(handles, labels, **legend_kwargs)


def rectangle(ax, x_range, y_range, **kwargs):
    """
    Draw rectangle on axis at specified ranges
    Args:
        ax:
        x_range:
        y_range:
        **kwargs:
    """
    ax.add_patch(
        patches.Rectangle((x_range[0], y_range[0]), x_range[1] - x_range[0],
                          y_range[1] - y_range[0], **kwargs)
    )


def get_2d_roc(x, y, matches, n_div=500, match_style='either'):
    """
    Get ROC curve for 2D data with two independent thresholds.
    Args:
        x: array-like[float]
            Length num_evidence array of evidence strength in x-dimension.
        y: array-like[float]
            Length num_evidence array of evidence strength in y-dimension.
        matches: array-like[bool]
            Length num_evidence array of bool, whether each point should be
            classified as True.
        n_div: int (Default=500)
            Number of percentiles (in x and y) to use when calculating the ROC
            curve.
        match_style: str ('either' or 'both')
            if 'either': a point is judged true if x > x_threshold OR
                         y > y_threshold.
            if 'both': a point is judged true if x > x_threshold AND
                       y > y_threshold.
    Returns:
        false_positives: numpy.array
            Array of false_positives along ROC curve, obtained by varying 2D
            thresholds.
        true_positives: numpy.array
            Array of true_positives along ROC curve, obtained by varying 2D
            thresholds.
        area_under_curve: float
            Area under ROC curve
    """
    assert match_style in ['either',
                           'both'], 'match_style must equal "either" or "both"'
    threshold_x = numpy.percentile(x, [q * 100.0 / (n_div - 1)
                                       for q in range(n_div)])
    threshold_y = numpy.percentile(y, [q * 100.0 / (n_div - 1)
                                       for q in range(n_div)])
    threshold_x = numpy.concatenate((numpy.array([float('-inf')]), threshold_x),
                                    axis=0)
    threshold_y = numpy.concatenate((numpy.array([float('-inf')]), threshold_y),
                                    axis=0)
    true_positives_2d = numpy.zeros((n_div, n_div), dtype=int)
    false_positives_2d = numpy.zeros((n_div, n_div), dtype=int)
    if match_style.lower() == 'either':
        for x, y, match in tqdm(zip(x, y, matches),
                                desc='finding ROC thresholds',
                                smoothing=0, mininterval=0.5, total=len(x)):
            i_x = bisect.bisect_right(threshold_x, x)
            i_y = bisect.bisect_right(threshold_y, y)
            if match:
                true_positives_2d[:i_x, :] += 1
                true_positives_2d[i_x:, :i_y + 1] += 1
            else:
                false_positives_2d[:i_x, :] += 1
                false_positives_2d[i_x:, :i_y + 1] += 1
    else:
        assert match_style.lower() == 'both', \
            'Unknown match_style (%s), must be "either" or "both"' % match_style
        for x, y, match in tqdm(zip(x, y, matches),
                                desc='finding ROC thresholds',
                                smoothing=0, mininterval=0.5, total=len(x)):
            i_x = bisect.bisect_right(threshold_x, x)
            i_y = bisect.bisect_right(threshold_y, y)
            if match:
                true_positives_2d[:i_x + 1, :i_y + 1] += 1
            else:
                false_positives_2d[:i_x, :i_y + 1] += 1

    false_positives_2d = false_positives_2d.reshape((false_positives_2d.size,)) / sum(~matches)
    true_positives_2d = true_positives_2d.reshape((true_positives_2d.size,)) / sum(matches)

    last_true_positives = float('-inf')
    last_false_positives = float('-inf')
    false_positives = []
    true_positives = []
    for x, y in sorted(zip(false_positives_2d, true_positives_2d),
                       key=operator.itemgetter(0)):
        if x <= last_false_positives:
            if y > last_true_positives:
                true_positives[-1] = y
                last_true_positives = y
        else:
            if y >= last_true_positives:
                false_positives.append(x)
                true_positives.append(y)
                last_true_positives = y
                last_false_positives = x
    false_positives = numpy.array(false_positives)
    true_positives = numpy.array(true_positives)

    area_under_curve = scipy.integrate.simps(true_positives, false_positives)
    return false_positives, true_positives, area_under_curve


def plot_interval_weights(_interval_weights, _coherent_weights, interval_matches,
                          weight_type, alpha=0.1):
    # make scatter plot
    fig1, ax1 = pyplot.subplots(1, 1)
    ax1.plot(_interval_weights[interval_matches],
             _coherent_weights[interval_matches],
             'cs', alpha=alpha,
             label='Good Intervals (n=%d)' % sum(interval_matches))
    ax1.plot(_interval_weights[~interval_matches],
             _coherent_weights[~interval_matches],
             'ro', alpha=alpha,
             label='Bad Intervals (n=%d)' % sum(~interval_matches))
    ax1.set_xlabel('Total Weight')
    ax1.set_ylabel('Coherent Weight')
    ax1.set_title('%s Weights' % weight_type)
    ax1.legend()


def add_colorbar(mappable, cax=None):
    ax = mappable.axes
    fig = ax.figure
    if cax is None:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def same_color_range(*mappables):
    vmin = min(m.get_clim()[0] for m in mappables)
    vmax = max(m.get_clim()[1] for m in mappables)
    for m in mappables:
        m.set_clim(vmin, vmax)


def hist1d(ax, x, bins=Default.num_1d_bins, xlabel="x", ylabel="counts", xscale=None, yscale="log", title=None,
           tick_labelrotation=Default.tick_labelrotation,
           x_is_good=None, x_max_left=numpy.inf, x_min_right=-numpy.inf, min_width=Default.min_1d_hist_width,
           histtype="stepfilled", orientation="vertical",
           title_fontsize=Default.title_fontsize, label_fontsize=Default.label_fontsize, **kwargs):
    if orientation == "horizontal":
        xlabel, ylabel = ylabel, xlabel
        xscale, yscale = yscale, xscale
    if isinstance(bins, int):
        num_bins = bins
        x_min = min(x.min(), x_max_left)
        x_max = max(x.max(), x_min_right)
        if x_max - x_min < min_width:
            x_min -= min_width / 2
            x_max += min_width / 2
        bins = numpy.linspace(x_min, x_max, num=num_bins + 1)
    # bins = "auto"
    if x_is_good is None:
        good_x = x
        bad_x = None
    else:
        good_x = x.compress(x_is_good)
        bad_x = x.compress(numpy.logical_not(x_is_good))

    bin_counts, bins, hist_patches = ax.hist(
        good_x, bins=bins, histtype=histtype, edgecolor='g', facecolor='g', orientation=orientation, **kwargs
    )
    if bad_x is not None:
        ax.hist(bad_x, bins=bins, histtype=histtype, edgecolor='r', facecolor='r', alpha=0.75,
                orientation=orientation, **kwargs)

    ax.tick_params(axis='x', labelrotation=tick_labelrotation)
    ax.autoscale(enable=True, tight=True)
    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)
    if title is not None:
        ax.set_title(title, fontsize=title_fontsize)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)

    return bin_counts, bins, hist_patches


def hist2d(x, y, x_bins=None, y_bins=None,
           xscale=None, yscale=None, zscale=None,
           xlabel: Optional[str] = 'x', ylabel: Optional[str] = 'y', zlabel='Counts',
           title='2D Histogram', mask=None, ax=None, cax=None,
           colorbar=True, cmap="viridis", num_2d_bins=Default.num_2d_bins,
           title_fontsize=Default.title_fontsize, label_fontsize=Default.label_fontsize):
    if x_bins is None:
        x_bins = min(max(10, int(len(x)**0.5)), num_2d_bins)
    if y_bins is None:
        y_bins = min(max(10, int(len(y)**0.5)), num_2d_bins)
    hist_mat, x_edges, y_edges = numpy.histogram2d(x, y, bins=[x_bins, y_bins])
    # rotate and flip hist_mat
    hist_mat = numpy.flipud(numpy.rot90(hist_mat))
    if mask is not None:
        # mark masked value as masked, they will not be plotted
        hist_mat = numpy.ma.masked_where(hist_mat == mask, hist_mat)
    if ax is None:
        fig, ax = pyplot.subplots(1, 1, num=next_fig())
    if zscale == "log":
        norm = colors.LogNorm()
    else:
        norm = None
    cmap = getattr(pyplot.cm, cmap)
    cmap.set_bad('k')
    mesh = ax.pcolormesh(x_edges, y_edges, hist_mat, norm=norm, cmap=cmap)
    xlim = (x_edges[0], x_edges[-1])
    ylim = (y_edges[0], y_edges[-1])
    if xscale is not None:
        ax.set_xscale(xscale)
        if xscale == "log":
            xlim = (min(x_edges[x_edges > 0]), max(x_edges[x_edges > 0]))
    if yscale is not None:
        ax.set_yscale(yscale)
        if yscale == "log":
            ylim = (min(y_edges[y_edges > 0]), max(y_edges[y_edges > 0]))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    if xlabel is not None and len(xlabel):
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if ylabel is not None and len(ylabel):
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if title is not None and len(title):
        ax.set_title(title, fontsize=title_fontsize)
    if colorbar:
        cbar = add_colorbar(mesh, cax=cax)
        # cbar = ax.figure.colorbar(mesh, ax=ax)
        cbar.set_label(zlabel)
    return mesh


def _sync_x(self, event):
    self.set_xlim(event.get_xlim(), emit=False)


def _sync_y(self, event):
    self.set_ylim(event.get_ylim(), emit=False)


# noinspection PyUnresolvedReferences
def sync_axes(axes: Sequence[pyplot.Axes], sync_x: bool = False, sync_y: bool = False, zoom_to_fit_all: bool = False):
    if not (sync_x or sync_y):
        return

    ax0 = axes[0]
    if zoom_to_fit_all:
        if sync_x:
            xlim = ax0.get_xlim()
            for ax2 in axes[1:]:
                xlim2 = ax2.get_xlim()
                xlim = (min(xlim[0], xlim2[0]), max(xlim[1], xlim2[1]))
            for ax in axes:
                ax.set_xlim(*xlim)
        if sync_y:
            ylim = ax0.get_ylim()
            for ax2 in axes[1:]:
                ylim2 = ax2.get_ylim()
                ylim = (min(ylim[0], ylim2[0]), max(ylim[1], ylim2[1]))
            for ax in axes:
                ax.set_ylim(*ylim)

    for ax in axes:
        if sync_x:
            ax.update_xlim = MethodType(_sync_x, ax)
        if sync_y:
            ax.update_ylim = MethodType(_sync_x, ax)

    for ax2 in axes[1:]:
        if sync_x:
            ax0.callbacks.connect("xlim_changed", ax2.update_xlim)
            ax2.callbacks.connect("xlim_changed", ax0.update_xlim)
        if sync_y:
            ax0.callbacks.connect("ylim_changed", ax2.update_ylim)
            ax2.callbacks.connect("ylim_changed", ax0.update_ylim)


def hist2d_fig(x, y, num_1d_bins=Default.num_1d_bins, num_2d_bins=Default.num_2d_bins,
               x_is_good=None, y_is_good=None,
               xlabel="x", ylabel="y", zlabel="z", title=None, tick_labelrotation=Default.tick_labelrotation,
               title_fontsize=Default.title_fontsize, label_fontsize=Default.label_fontsize,
               **kwargs):
    """construct a figure where a 2D histogram is flanked by 1D histograms for each axis"""
    fig = pyplot.figure(num=next_fig(), constrained_layout=True)
    gridspec = fig.add_gridspec(2, 3, height_ratios=[0.8, 0.2], width_ratios=[0.2, 0.75, 0.05])

    # make top-right 2d histogram: z(x, y)
    ax0 = fig.add_subplot(gridspec[0, 1])
    cax = fig.add_subplot(gridspec[0, 2])
    if title is None:
        title = "%s (%s, %s)" % (zlabel, xlabel, ylabel)
    hist2d(
        x, y, ax=ax0, xlabel=None, ylabel=None, zlabel=zlabel,
        title=title, num_2d_bins=num_2d_bins, cax=cax, **kwargs
    )
    ax0.set_xticklabels([])
    ax0.set_yticklabels([])

    hist1d_kwargs = {k: v for k, v in kwargs.items() if not k.startswith("z")}
    # make left 1d histogram: z(y)
    ax1 = fig.add_subplot(gridspec[0, 0])
    hist1d(
        ax1, y, x_is_good=y_is_good, bins=num_1d_bins, xlabel=ylabel, ylabel=zlabel, title=None,
        tick_labelrotation=tick_labelrotation, orientation="horizontal", **hist1d_kwargs,
        title_fontsize=title_fontsize, label_fontsize=label_fontsize
    )

    # make bottom 1d histogram: z(x)
    ax2 = fig.add_subplot(gridspec[1, 1])
    hist1d(
        ax2, x, x_is_good=x_is_good, bins=num_1d_bins, xlabel=xlabel, ylabel=zlabel, title=None,
        tick_labelrotation=tick_labelrotation, title_fontsize=title_fontsize, label_fontsize=label_fontsize,
        **hist1d_kwargs
    )

    sync_axes([ax0, ax1], sync_x=False, sync_y=True)
    sync_axes([ax0, ax2], sync_x=True, sync_y=False)

    # ax0.update_xlim = MethodType(_sync_x, ax0)
    # ax0.update_ylim = MethodType(_sync_y, ax0)
    # ax1.update_ylim = MethodType(_sync_y, ax1)
    # ax2.update_xlim = MethodType(_sync_x, ax2)
    #
    # ax0.callbacks.connect("ylim_changed", ax1.update_ylim)
    # ax0.callbacks.connect("xlim_changed", ax2.update_xlim)
    # ax1.callbacks.connect("ylim_changed", ax0.update_ylim)
    # ax2.callbacks.connect("xlim_changed", ax0.update_xlim)


def save_figures(figure_save_file: str, *figures: Optional[pyplot.Figure]):
    """
    Save one or more figures into a single pdf
    Args:
        figure_save_file: str
            Path to save pdf
        *figures: pyplot.Figure or None
            Additional arguments are saved as figures in the pdf.
            As a convenience, this function skips figures that are None
    """
    with PdfPages(figure_save_file) as pdf:
        for fig in figures:
            if fig is None:
                continue
            pdf.savefig(fig, bbox_inches="tight")
