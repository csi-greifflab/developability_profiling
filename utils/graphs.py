"""
Functions to generate graphs
version 1.0

Author: Matteo Pariset

"""

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

# plt.style.use('ggplot')
# plt.rcParams["font.family"] = "Helvetica"

def to_color(seq):
    seq -= seq.min()
    if np.isclose(seq.max(), 0):
        return np.ones_like(seq)
    else:
        seq /= seq.max()
        return seq

def to_centered_unit_interval(seq, median, vmin=None, vmax=None):
    seq -= median

    high_samples = seq > 0
    low_samples = seq <= 0

    if high_samples.sum() != 0:
        if vmax is None:
            high_width = seq[high_samples].max()
        else:
            high_width = vmax - median

        if not np.isclose(high_width, 0):
            seq[high_samples] /= (2*high_width)

    if low_samples.sum() != 0:
        if vmin is None:
            low_width = -seq[low_samples].min()
        else:
            low_width = median - vmin

        if not np.isclose(low_width, 0):
            seq[low_samples] /= (2*low_width)

    return seq

def to_grey_zone_mask(seq, grey_zone):
    return np.logical_and(seq > -grey_zone/2, seq < grey_zone/2)

class MedianCenteredNorm(mpl.colors.Normalize):
    def __init__(self, median, vmin, vmax, grey_zone, clip=False) -> None:
        self.median = median
        self.grey_zone = grey_zone
        super().__init__(vmin, vmax, clip=clip)

    def __call__(self, value, clip=None):
        if isinstance(value, float):
            value = np.array([value])
        # else:
        #     value = np.copy(value)

        grey_zone_low = self.median - (self.median - self.vmin) * self.grey_zone / 2
        grey_zone_high = self.median + (self.vmax - self.median) * self.grey_zone / 2

        x = [self.vmin, grey_zone_low, grey_zone_high, self.vmax]
        y = [0, .5 - self.grey_zone/2, .5 + self.grey_zone/2, 1]

        return np.interp(value, x, y)

    def inverse(self, value):
        return super().inverse(value)


def to_colormap(seq, median, grey_zone=0., colormap='cividis'):
    norm = MedianCenteredNorm(median=median, vmin=seq.min(), vmax=seq.max(), grey_zone=grey_zone)

    low_cm = mpl.cm.get_cmap(colormap)(np.linspace(0, .5 - grey_zone/2, int((1-grey_zone)*256)))
    grey_slice_ptxs = int(grey_zone*512)
    grey_cm = np.hstack((np.zeros((grey_slice_ptxs, 3)), np.ones((grey_slice_ptxs, 1))))
    high_cm = mpl.cm.get_cmap(colormap)(np.linspace(.5 + grey_zone/2, 1, int((1-grey_zone)*256)))

    blacked_out_median_cm = mpl.colors.LinearSegmentedColormap.from_list('blacked_out_median', np.vstack((low_cm, grey_cm, high_cm)))

    sc_mapper = mpl.cm.ScalarMappable(norm=norm, cmap=blacked_out_median_cm)

    seq = to_centered_unit_interval(seq, median)

    # Black-out the median
    seq[to_grey_zone_mask(seq, grey_zone)] = np.nan

    # Bring to [0,1]
    seq += .5

    cmap = mpl.cm.get_cmap(colormap)
    return np.array([cmap(x) for x in seq]), sc_mapper

def as_title(dim_reduction, embedding_model, pca_computation, pca_evaluation, embs_reduction, label_type, first_axis, second_axis, line_sep="\n"):
    title = f"{dim_reduction.__class__.__name__.lower()}-{embedding_model}-{embs_reduction}{line_sep}{pca_computation}-{first_axis},{second_axis}{line_sep}{pca_evaluation}-{label_type}"
    return title

def _color_based_on_median(labels, min_alpha=0., max_alpha=1., median=None, **kwargs):
    cs = np.zeros((labels.shape[0], 4))
    if median is None:
        median = labels.median()
    cs[:,:], cs_mapper = to_colormap(labels, median, **kwargs)
    cs[:,3] = np.clip(cs[:,3], min_alpha, max_alpha)
    return cs, cs_mapper

# Double-shuffle-with-smudge coloring depth computation
def compute_depth_ordering(labels, grey_zone=0., median=None, **kwargs):
    if labels.dtype.name == "category":
        return np.random.permutation(labels.shape[0])

    if median is None:
        median = labels.median()

    seq = to_centered_unit_interval(labels, median)

    gz_mask = to_grey_zone_mask(seq.values, grey_zone)
    gz_ordering = np.argsort(gz_mask*1)[::-1]
    
    boundary_idx = gz_mask.sum()
    smudge_boundary_idx = np.maximum(int(boundary_idx * .7), 0)

    # Sanity checks
    assert boundary_idx == 0 or gz_mask[gz_ordering[:boundary_idx]].mean() == 1., f"Found 0s among extreme labels, {gz_mask[gz_ordering[:boundary_idx]].mean()}"
    assert gz_mask[gz_ordering[boundary_idx:]].mean() == 0., f"Found 1s among non-extreme labels, {gz_mask[gz_ordering[boundary_idx:]].mean()}"
    
    gz_ordering[:boundary_idx] = np.random.permutation(gz_ordering[:boundary_idx])
    gz_ordering[smudge_boundary_idx:] = np.random.permutation(gz_ordering[smudge_boundary_idx:])

    return gz_ordering

def _color_based_on_category(labels, alpha=.65, na_alpha=.05, category_cmap='tab20', **kwargs):
    cs = np.zeros((labels.shape[0], 4))

    cs[:,:] = mpl.cm.get_cmap(category_cmap)(labels.cat.codes)
    cs[:,3] = alpha

    # Handle NA
    na_idxs = labels.isna()
    cs[na_idxs,:3] = 0.
    cs[na_idxs,3] = na_alpha

    return cs

def get_colors_from_labels(labels, *args, **kwargs):
    # Check if labels are categorical
    if labels.dtype.name == "category":
        return _color_based_on_category(labels, *args, **kwargs), None
    else:
        return _color_based_on_median(labels, *args, **kwargs)

def show_legend(c_plt):
    leg = c_plt.legend()
    for i in range(len(leg.legendHandles)):
        # Good hack from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
        leg.legendHandles[i]._sizes = [10]

    return leg

def plot_pca(embs, labels, first_axis, second_axis, invisible_zone=0., grey_zone=0., labels_preview=False, mask=None, ax=(plt, plt), c=None, label=None, colormap='RdYlGn', min_alpha=0., max_alpha=1., **kwargs):
    if np.isclose(invisible_zone, 0.):
        extreme_idxs = np.ones(labels.shape[0]).astype(bool)
    else:
        extreme_idxs = np.logical_not(((labels < labels.median() + invisible_zone*(labels.max()-labels.median())) * (labels > labels.median() - invisible_zone*(labels.median()-labels.min()))).to_numpy())

    if mask is not None:
        extreme_idxs = np.logical_and(extreme_idxs, mask)

    extreme_labels = labels[extreme_idxs]

    if labels_preview:
        ax[1].scatter(np.arange(extreme_labels.shape[0]), extreme_labels)
        ax[0].show()

    cbar = None
    
    if grey_zone < 1.:
        if c is None:
            cs, cs_mapper = get_colors_from_labels(extreme_labels, grey_zone=grey_zone, min_alpha=min_alpha, max_alpha=max_alpha, colormap=colormap)
            cs = cs[extreme_idxs]

            if cs_mapper is not None:
                labels_median = extreme_labels.median()
                cbar = ax[0].colorbar(cs_mapper, spacing='proportional')
                cbar.set_ticks(cbar.get_ticks()[1:-1])
        else:
            cs = c
        depth_order = compute_depth_ordering(extreme_labels, grey_zone=grey_zone)

        if labels.dtype.name == "category":
            unique_cats = extreme_labels.reset_index(drop=True).reset_index(drop=False).drop_duplicates(subset=extreme_labels.name)
            unique_cats.join(unique_cats[extreme_labels.name].cat.codes.rename(extreme_labels.name), rsuffix="_code").reset_index().sort_values(extreme_labels.name).apply(
                lambda cat: ax[1].scatter(
                    embs[extreme_idxs, first_axis][int(cat['index'])],
                    embs[extreme_idxs, second_axis][int(cat['index'])], 
                    alpha=max_alpha, 
                    label=cat[extreme_labels.name],
                    c=[mpl.cm.get_cmap('tab20')(cat[extreme_labels.name + "_code"])[:-1] if cat[extreme_labels.name + "_code"] >= 0 else (0, 0, 0, .1)],
                    s=1,
                    zorder=0), 
                axis=1)
            show_legend(ax[1])
        ax[1].scatter(embs[extreme_idxs,first_axis][depth_order], embs[extreme_idxs, second_axis][depth_order], c=cs[depth_order], **kwargs)
    else:
        ax[1].scatter(embs[extreme_idxs,first_axis], embs[extreme_idxs, second_axis], c=[[.85, .85, .85, .6]], **kwargs)

    if label is not None:
        for idx, t in enumerate(label):
            ax[1].annotate(t, (embs[idx, first_axis], embs[idx, second_axis]))

    return cbar
