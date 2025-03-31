import nibabel as nib
from gradec.decode import LDADecoder, TermDecoder
from gradec.fetcher import _fetch_classification, _fetch_features, _fetch_frequencies
from gradec.plot import plot_cloud, plot_radar
from gradec.utils import _decoding_filter



# get st_cen
cen_idx = 2
cen_file = r"C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r400_ro1\regional_bag_st_cen_2in3.dscalar.nii"

cen_cii = nib.load(cen_file)
cen_map = cen_cii.get_fdata()


###############################################################################
# Decoding
# -----------------------------------------------------------------------------
out_dir = r"C:\Users\lihon\Desktop\abcd_study\pnc\meta_analysis"
data_dir = r"C:\Users\lihon\Desktop\abcd_study\pnc\meta_analysis\gradec"

SPACE, DENSITY = "fsLR", "32k"
DSET = "neuroquery"
MODEL = "term"

# Fit the decoder
calc_p = False

decode = TermDecoder(space=SPACE, density=DENSITY, calc_pvals=calc_p, data_dir=data_dir)
decode.fit(DSET)

# Load features for visualization
features = _fetch_features(DSET, MODEL, data_dir=data_dir)
classification, class_lst = _fetch_classification(DSET, MODEL, data_dir=data_dir)

frequencies = None
class_lst = None

corrs_df = decode.transform(cen_map, method="correlation")
if calc_p:
    corrs_df, pvals_df, corr_pvals_df = corrs_df


###############################################################################
# Visualization
# -----------------------------------------------------------------------------
pos_cor = True

for map_i, (c_map, segment) in enumerate(zip(cen_map, corrs_df.columns)):
    seg_df = corrs_df[[segment]]
    # Filter results
    if not pos_cor:
        seg_df *= -1
    
    filtered_lst = _decoding_filter(
        seg_df,
        features,
        classification,
        pos_corr=pos_cor,
        freq_by_topic=frequencies,
        class_by_topic=class_lst,
        class_to_keep=["Functional"]
    )
    filtered_lst = list(filtered_lst)
    if len(filtered_lst)==2:
        filtered_lst.append(None)

    filtered_df, filtered_features, filtered_frequencies = filtered_lst
    
    filtered_df.columns = ["r"]
    corrs = filtered_df["r"].to_numpy()

    # Word cloud plot
    cloud_fig = plot_cloud(
        corrs,
        filtered_features,
        MODEL,
        frequencies=filtered_frequencies,
        dpi=300
    )
    out_fig = out_dir + r"\wordcloud_" + DSET + "_" + MODEL + "_" + str(pos_cor) + ".jpg"
    cloud_fig.savefig(out_fig, bbox_inches="tight", dpi=300)
    
    # radar plot
    radar_fig = plot_radar(
        corrs,
        filtered_features,
        MODEL
    )
    out_fig = out_dir + r"\radar_" + DSET + "_" + MODEL + "_" + str(pos_cor) + ".jpg"
    radar_fig.savefig(out_fig, bbox_inches="tight", dpi=300)

print("Done!")


