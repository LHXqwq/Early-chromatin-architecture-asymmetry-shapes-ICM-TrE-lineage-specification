import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["USE_TF"] = "0"

import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
import hdbscan
from umap import UMAP

# 设置绘图样式
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['svg.hashsalt'] = 'hello'

# =============================
# 读取数据
# =============================
data = pd.read_csv("IVF_350_Ps_curve_new.txt", sep="\t")
data['ID'] = data['ID'].str.strip()

# =============================
# 横坐标
# =============================
x = data.columns[1:].astype(float)
x_new = np.linspace(x.min(), x.max(), 500)

# =============================
# feature extraction
# =============================
def extract_feature(row):
    y = row.values[1:].astype(float)

    if not np.isfinite(y).all():
        return None

    try:
        spl = make_interp_spline(x, y, k=3)
        y_smooth = spl(x_new)

        mask = (x_new >= 5e5) & (x_new <= 2e7)
        y_use = y_smooth[mask]

        y_use = np.clip(y_use, 1e-10, None)
        y_log = np.log(y_use)

        y_log = (y_log - np.mean(y_log)) / (np.std(y_log) + 1e-8)

        return y_log

    except:
        return None


# =============================
# UMAP + plotting style (3D version)
# =============================
def plot_umap_3d(features, labels, title, save_path=None):

    reducer = UMAP(
        n_components=3,
        n_neighbors=15,
        min_dist=0.2,
        metric='euclidean',
        random_state=42
    )

    embedding = reducer.fit_transform(features)
    labels = np.array(labels)

    color_map = {
        1: "#3cb44b",
        2: "#e6194b",
        0: "#e6194b",
        -1: "gray"
    }

    plt.rcParams.update({
        "font.size": 22,
        "axes.titlesize": 32,
        "axes.labelsize": 22,
        "legend.fontsize": 20
    })

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    for lab in np.unique(labels):
        idx = labels == lab

        ax.scatter(
            embedding[idx, 0],
            embedding[idx, 1],
            embedding[idx, 2],
            s=70,
            alpha=0.85,
            c=color_map.get(lab, "black"),
            label=str(lab)
        )

    ax.view_init(elev=20, azim=135)

    ax.set_title(title, pad=30)
    ax.set_xlabel("UMAP 1", labelpad=20)
    ax.set_ylabel("UMAP 2", labelpad=20)
    ax.set_zlabel("")

    ax.text2D(
        1.02, 0.5, "UMAP 3",
        transform=ax.transAxes,
        fontsize=20,
        rotation=90,
        va='center'
    )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.legend(loc='upper right', frameon=False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)

    plt.show()


# =============================
# Hierarchical clustering
# =============================
def run_hierarchical(stage_prefix, stage_name):

    subset = data[data['ID'].str.startswith(stage_prefix)]

    features, cell_ids = [], []

    for _, row in subset.iterrows():
        feat = extract_feature(row)
        if feat is not None:
            features.append(feat)
            cell_ids.append(row['ID'])

    features = np.array(features)

    Z = linkage(features, method='ward')
    labels = fcluster(Z, t=2, criterion='maxclust')

    pd.DataFrame({
        "ID": cell_ids,
        "Cluster": labels
    }).to_csv(f"{stage_name}_hierarchical_ps.txt", sep="\t", index=False)

    plot_umap_3d(features, labels, f"{stage_name} Hierarchical + UMAP (ps curve)",
                 f"{stage_name}_hierarchical_umap_ps.svg")


# =============================
# HDBSCAN clustering
# =============================
import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
import hdbscan
import matplotlib.pyplot as plt
import umap
import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["USE_TF"] = "0"

# =============================
# 读取数据
# =============================
data = pd.read_csv("IVF_350_Ps_curve_new.txt", sep="\t")
data['ID'] = data['ID'].str.strip()

x = data.columns[1:].astype(float)
x_new = np.linspace(x.min(), x.max(), 500)

# =============================
# feature extraction
# =============================
def extract_feature(row):

    y = row.values[1:].astype(float)

    if not np.isfinite(y).all():
        return None

    try:
        spl = make_interp_spline(x, y, k=3)
        y_smooth = spl(x_new)

        mask = (x_new >= 5e5) & (x_new <= 2e7)
        y_use = y_smooth[mask]

        y_use = np.clip(y_use, 1e-10, None)
        y_log = np.log(y_use)

        y_log = (y_log - np.mean(y_log)) / (np.std(y_log) + 1e-8)

        return y_log

    except:
        return None


# =============================
# HDBSCAN → 强制2分类函数
# =============================
def run_hdbscan_force2(stage_prefix, stage_name):

    subset = data[data['ID'].str.startswith(stage_prefix)]

    features, cell_ids = [], []

    for _, row in subset.iterrows():
        feat = extract_feature(row)
        if feat is not None:
            features.append(feat)
            cell_ids.append(row['ID'])

    features = np.array(features)

    # =============================
    # step 1: PCA
    # =============================
    if features.shape[0] > 10:
        features = PCA(n_components=10).fit_transform(features)

    # =============================
    # step 2: HDBSCAN (over-cluster)
    # =============================
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=5,
        min_samples=1,
        cluster_selection_method='leaf'
    )

    labels = clusterer.fit_predict(features)

    # =============================
    # step 3: 去掉 noise
    # =============================
    mask = labels != -1

    if mask.sum() < 10:
        print(f"{stage_name}: too few non-noise points")
        return

    X = features[mask]
    lab = labels[mask]

    unique_clusters = np.unique(lab)

    # =============================
    # step 4: cluster centroid
    # =============================
    centroids = []
    for c in unique_clusters:
        centroids.append(X[lab == c].mean(axis=0))

    centroids = np.array(centroids)

    # =============================
    # step 5: 强制2分类
    # =============================
    agg = AgglomerativeClustering(n_clusters=2)
    coarse = agg.fit_predict(centroids)

    cluster_map = {old: new for old, new in zip(unique_clusters, coarse)}

    final_labels = np.array([cluster_map[i] for i in lab])

    # =============================
    # save result
    # =============================
    out_df = pd.DataFrame({
        "ID": np.array(cell_ids)[mask],
        "Cluster": final_labels
    })

    out_df.to_csv(f"{stage_name}_hdbscan_force2_ps.txt", sep="\t", index=False)

    # =============================
    # UMAP
    # =============================
    plot_umap_3d(
        X,
        final_labels,
        f"{stage_name} HDBSCAN + UMAP (ps curve)",
        f"{stage_name}_hdbscan_umap_ps.svg"
    )


# =============================
# RUN BOTH STAGES
# =============================
run_hierarchical("IVF-2-", "IVF-2C")
run_hierarchical("IVF-4-", "IVF-4C")

run_hdbscan_force2("IVF-2-", "IVF-2C")
run_hdbscan_force2("IVF-4-", "IVF-4C")