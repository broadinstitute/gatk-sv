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
from typing import Optional, Union, Tuple, List, Sequence, Mapping, Text, Any, Callable
from types import MethodType, MappingProxyType


Numeric = Union[int, float, numpy.integer, numpy.floating]
Record = Union[numpy.record, Tuple]
Vector = Union[pandas.Series, numpy.ndarray]
ListLike = Union[List, pandas.Series, numpy.ndarray]


class Default:
    tick_labelrotation = -30
    num_1d_bins = 400
    num_2d_bins = 200
    min_1d_hist_width = 1.0e-6
    title_fontsize = 16
    label_fontsize = 16
    tick_fontsize = 12
    num_palette_colors = 12
    color_palette = "colorblind"
    # Default params overrides for matplotlib
    matplotlib_params = MappingProxyType({
        'savefig.dpi': 150,  # to adjust notebook inline plot size
        'axes.labelsize': label_fontsize,
        'axes.titlesize': title_fontsize,
        'font.size': label_fontsize,
        'legend.fontsize': label_fontsize,
        'xtick.labelsize': tick_fontsize,
        'ytick.labelsize': tick_fontsize,
        'figure.figsize': [9.75, 6]
    })
    use_seaborn = True
    mask_color = 'k'
    hist2d_colormap = "viridis"


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


def ordered_legend(ax: pyplot.Axes, labels_in_order: Sequence, **legend_kwargs) -> matplotlib.legend.Legend:
    """ get legend, specifying the order of the labels """
    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: labels_in_order.index(t[0])))
    return ax.legend(handles, labels, **legend_kwargs)


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


def barh(ax, x_baseline, edges, values, **kwargs):
    y_low = edges[0]
    for value, y_high in zip(values, edges[1:]):
        # rectangle(ax, [x_baseline, x_baseline + value], [edges[0], edges[1]])
        if numpy.isfinite(value):
            ax.add_patch(
                patches.Rectangle(
                    (x_baseline, y_low), value, y_high - y_low, **kwargs
                )
            )
        y_low = y_high


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
        false_positives: numpy.ndarray
            Array of false_positives along ROC curve, obtained by varying 2D
            thresholds.
        true_positives: numpy.ndarray
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


def _fast_hist2d(
        x: numpy.ndarray,
        y: numpy.ndarray,
        x_edges: numpy.ndarray,
        y_edges: numpy.ndarray
):
    num_x_bins = len(x_edges) - 1
    num_y_bins = len(y_edges) - 1
    max_ravel_bin = num_x_bins * num_y_bins - 1
    bin_type = numpy.min_scalar_type(max_ravel_bin)
    eps = 0 if numpy.issubdtype(x.dtype, numpy.integer) \
        else numpy.finfo(x.dtype).resolution

    def _avg_diff(_edges) -> float:
        return (_edges[-1] - _edges[0]) / (len(_edges) - 1)

    def _is_uniform_spaced(_edges: numpy.ndarray) -> bool:
        return numpy.allclose(_avg_diff(_edges), numpy.diff(_edges))

    def _get_bin_indices(_values: numpy.ndarray, _edges: numpy.ndarray) -> numpy.ndarray:
        """ Get indices corresponding to regularly-spaced bins """
        # pad bins with eps to avoid rounding errors leading to invalid bin assignments (<=0 or >=num_bins)
        _delta = 2 * eps + _avg_diff(_edges)
        _base = _edges[0] - eps
        return numpy.floor(
            (_values - _base) / _delta
        ).astype(bin_type)

    def _get_ravel_indices():
        # want to flip y and rotate 90 at *THIS* stage so we don't have to do that operation to hist_mat
        return _get_bin_indices(x, x_edges) + num_x_bins * _get_bin_indices(y, y_edges)

    if _is_uniform_spaced(x_edges) and _is_uniform_spaced(y_edges):
        # can do faster 2d histogram
        ravel_indices = _get_ravel_indices()
        bin_counts = numpy.bincount(ravel_indices)
        hist_mat = numpy.zeros((num_x_bins, num_y_bins), dtype=bin_type)
        hist_mat.ravel()[:len(bin_counts)] = bin_counts
    else:
        # have to do it slowly
        hist_mat, x_edges, y_edges = numpy.histogram2d(x, y, bins=[x_edges, y_edges], density=False)
        # rotate and flip hist_mat
        hist_mat = numpy.flipud(numpy.rot90(hist_mat))
    return hist_mat, x_edges, y_edges


def hist2d(
        ax: pyplot.Axes,
        x: Union[numpy.ndarray, pandas.Series],
        y: Union[numpy.ndarray, pandas.Series],
        x_bins: Union[Sequence[int], int, None] = None,
        y_bins: Union[Sequence[int], int, None] = None,
        xscale: Optional[str] = None,
        yscale: Optional[str] = None,
        zscale: Optional[str] = None,
        xlabel: Optional[str] = 'x',
        ylabel: Optional[str] = 'y',
        zlabel: Optional[str] = 'Counts',
        xlim: Optional[Sequence[float]] = None,
        ylim: Optional[Sequence[float]] = None,
        title='2D Histogram',
        mask_condition: Optional[Callable[[numpy.ndarray], numpy.ndarray]] = None,
        mask_color: str = Default.mask_color,
        colorbar_axes: Optional[pyplot.Axes] = None,
        colorbar: bool = True,
        cmap: Union[str, matplotlib.colors.Colormap] = Default.hist2d_colormap,
        num_2d_bins: int = Default.num_2d_bins,
        title_fontsize: int = Default.title_fontsize,
        label_fontsize: int = Default.label_fontsize,
        tick_fontsize: int = Default.tick_fontsize
) -> (Optional[matplotlib.collections.QuadMesh], numpy.ndarray, numpy.ndarray, numpy.ndarray):
    f"""
    Plot a 2D histogram (heat map) on an axes
    Args:
        ax: pyplot.Axes
            Axes to draw on
        x: numpy.ndarray or pandas.Series
            Array of x values
        y: numpy.ndarray, or pandas.Series
            Array of y values
        x_bins: array-like, int, or None = None
            if array, use directly as bin edges
            if an int, specify number of bins to use for x coordinate, and space bins based on range and xscale
            If None, try to select a reasonable number of bins automatically.
        y_bins:  array-like, int, or None = None
            if array, use directly as bin edges
            if an int, specify number of bins to use for y coordinate, and space bins based on range and yscale
            If None, try to select a reasonable number of bins automatically.
        xscale: Optional[str] = None
            String to pass through to ax.set_xscale(), e.g. "log", "symlog"
        yscale: Optional[str] = None
            String to pass through to ax.set_yscale(), e.g. "log", "symlog"
        zscale: Optional[str] = NOne
            Either "log" (log-scale the colors) or None.
        xlabel: Optional[str] = None
            label for x axis
        ylabel: Optional[str] = None
            label for y axis
        zlabel: Optional[str] = None
            label for intensity / colorbar
        xlim: Optional[Sequence[float]] = None
            If list/tuple of (min, max) values supplied, use these limits for x axis. Otherwise set from data range.
        ylim: Optional[Sequence[float]] = None
            If list/tuple of (min, max) values supplied, use these limits for x axis. Otherwise set from data range.
        title: Optional[str] = None
            title for axes
        mask_condition: Optional[Callable] = None
            If provided, should be a function that takes a matrix of histogram values as input, and outputs a matrix of
            bool that is True everywhere that should be masked. e.g. `mask_condition = lambda x: x == my_masked_value`
            will result in a plot that has values masked (set to mask_color) whenever the counts == my_masked_value.
        mask_color: str = {Default.mask_color}
            Color with which to paint masked regions of the heat map
        colorbar_axes: Optional[pyplot.Axes] = None
            If non-None, draw colorbar onto this axes
        colorbar: bool = True
            If True, draw a colorbar
        cmap: str = {Default.hist2d_colormap}
        num_2d_bins: int = {Default.num_2d_bins}
            Maximum number of bins for each axis, used when auto-determining number of bins
        title_fontsize: int = {Default.title_fontsize}
            Font size for title
        label_fontsize: int = {Default.label_fontsize}
            Font size for axis labels
        tick_fontsize: int = {Default.tick_fontsize}
    Returns:
        mesh: matplotlib.Collections.QuadMesh
            matplotlib object containing heatmap. Will be None if no drawable data is provided.
        hist_mat: numpy.ndarray
            num_x_bins x num_y_bins matrix of counts for each bin on the grid
        x_edges: numpy.ndarray
            array of length num_x_bins + 1, with the edges of each x bin
        y_edges: numpy.ndarray
            array of length num_y_bins + 1, with the edges of each y bin
    """
    if isinstance(x, pandas.Series):
        x = x.values
    if isinstance(y, pandas.Series):
        y = y.values
    assert len(x) == len(y)
    is_finite = numpy.logical_and(numpy.isfinite(x), numpy.isfinite(y))
    x = x.compress(is_finite)
    y = y.compress(is_finite)
    num_points = len(x)
    # it's faster to work with uniform bins, so transform data first, then scale edges for plot
    if xscale == "log":
        x = numpy.log10(x)
    if yscale == "log":
        y = numpy.log10(y)

    def _get_edges_arr(
            _bins: Union[None, int, Sequence[int]], _scale: Optional[str], _range: Tuple[float, float]
    ) -> numpy.ndarray:
        if _bins is None:
            # want num_x_bins * num_y_bins ~ max(100, sqrt(num_points))
            _bins = min(
                max(10, int(num_points ** 0.25)),
                num_2d_bins
            )
        # if _bins was supplied as a sequence, transform it into appropriate space, otherwise make uniform bins
        if isinstance(_bins, int):
            _bins = numpy.linspace(*_range, _bins + 1)
        elif _scale == "log":
            _bins = numpy.log10(_bins)
        else:
            _bins = numpy.array(_bins)
        return _bins

    if len(x) > 0:
        # we have some sensible data
        # noinspection PyTypeChecker
        x_edges = _get_edges_arr(x_bins, xscale, (numpy.nanmin(x), numpy.nanmax(x)))
        # noinspection PyTypeChecker
        y_edges = _get_edges_arr(y_bins, yscale, (numpy.nanmin(y), numpy.nanmax(y)))
        hist_mat, x_edges, y_edges = _fast_hist2d(x, y, x_edges=x_edges, y_edges=y_edges)
        # transform edge coordinates back to correct values
        if xscale == "log":
            x_edges = 10 ** x_edges
        if yscale == "log":
            y_edges = 10 ** y_edges

        ax.grid(False)
        mesh = ax.pcolormesh(
            x_edges, y_edges,
            hist_mat if mask_condition is None else numpy.ma.masked_where(mask_condition(hist_mat), hist_mat),
            norm=colors.LogNorm() if zscale == "log" else colors.Normalize(), cmap=cmap, edgecolors="face"
        )
    else:
        hist_mat, x_edges, y_edges = numpy.array([]).reshape((0, 0)), numpy.array([]), numpy.array([])
        mesh = None

    if not isinstance(cmap, matplotlib.colors.Colormap):
        cmap = getattr(pyplot.cm, cmap).copy()
    if mask_condition is not None:
        cmap.set_bad(mask_color)
    if xlim is None and len(x_edges) > 0:
        xlim = (min(x_edges[x_edges > 0]), max(x_edges[x_edges > 0])) if xscale == "log" else (x_edges[0], x_edges[-1])
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is None and len(y_edges) > 0:
        ylim = (min(y_edges[y_edges > 0]), max(y_edges[y_edges > 0])) if yscale == "log" else (y_edges[0], y_edges[-1])
    if ylim is not None:
        ax.set_ylim(*ylim)
    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)
    if xlabel is not None and len(xlabel):
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if ylabel is not None and len(ylabel):
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if title is not None and len(title):
        ax.set_title(title, fontsize=title_fontsize)
    if tick_fontsize is not None:
        ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
    if colorbar and mesh is not None:
        # noinspection PyArgumentList
        cbar = add_colorbar(mesh, ax=[ax]) if colorbar_axes is None \
            else add_colorbar(mesh, cax=colorbar_axes)
        if zlabel is not None:
            cbar.set_label(zlabel)
    return mesh, hist_mat, x_edges, y_edges


def _sync_x(self, event):
    """ Event handler for updates to x axis: set x axis to desired limits """
    self.set_xlim(event.get_xlim(), emit=False)


def _sync_y(self, event):
    """ Event handler for updates to y axis: set y axis to desired limits """
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
        title=title, num_2d_bins=num_2d_bins, colorbar_axes=cax, **kwargs
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
