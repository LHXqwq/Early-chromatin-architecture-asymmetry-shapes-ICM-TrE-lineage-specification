import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from umap import UMAP
import matplotlib.pyplot as plt

# 设置绘图样式
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['svg.hashsalt'] = 'hello'

# =========================
# 1. 读取数据
# =========================
df = pd.read_csv("IVF_scAB_eigen.tsv", sep="\t")
df["bin_id"] = df["chrom"].astype(str) + "_" + df["bin"].astype(str)

matrix = df.pivot_table(
    index="cell",
    columns="bin_id",
    values="scAB_eigen"
).fillna(0)

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



# =========================
# 3. Hierarchical clustering（强制2类）
# =========================
def run_hierarchical(stage_prefix, name):

    sub = matrix[matrix.index.str.startswith(stage_prefix)]
    X = sub.values

    model = AgglomerativeClustering(n_clusters=2)
    labels = model.fit_predict(X)

    out = pd.DataFrame({
        "cell": sub.index,
        "cluster": labels
    })

    out.to_csv(f"{name}_hierarchical_ab.tsv", sep="\t", index=False)

    plot_umap_3d(
        X,
        labels,
        f"{name} Hierarchical + UMAP (scAB)",
        f"{name}_hierarchical_umap_ab.svg"
    )


# =========================
# RUN
# =========================
run_hierarchical("IVF-2-", "IVF-2C")
run_hierarchical("IVF-4-", "IVF-4C")

import pandas as pd
import numpy as np
import hdbscan
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from umap import UMAP
import matplotlib.pyplot as plt

# =========================
# 1. 读取数据
# =========================
df = pd.read_csv("IVF_scAB_eigen.tsv", sep="\t")
df["bin_id"] = df["chrom"].astype(str) + "_" + df["bin"].astype(str)

matrix = df.pivot_table(
    index="cell",
    columns="bin_id",
    values="scAB_eigen"
).fillna(0)

# =========================
# 2. UMAP绘图（同风格）
# =========================
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


# =========================
# 3. HDBSCAN → 强制2类
# =========================
def run_hdbscan(stage_prefix, name):

    sub = matrix[matrix.index.str.startswith(stage_prefix)]
    X = sub.values

    # PCA（稳定）
    X = PCA(n_components=10).fit_transform(X)

    # HDBSCAN
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=5,
        min_samples=1,
        cluster_selection_method='leaf'
    )

    labels = clusterer.fit_predict(X)

    # 去noise
    mask = labels != -1
    X2 = X[mask]
    lab = labels[mask]

    if len(np.unique(lab)) < 2:
        print(f"{name}: not enough clusters")
        return

    # cluster → centroid → 强制2类
    centroids = []
    uniq = np.unique(lab)

    for c in uniq:
        centroids.append(X2[lab == c].mean(axis=0))

    centroids = np.array(centroids)

    agg = AgglomerativeClustering(n_clusters=2)
    coarse = agg.fit_predict(centroids)

    map_dict = {old: new for old, new in zip(uniq, coarse)}

    final = np.array([map_dict[i] for i in lab])

    out = pd.DataFrame({
        "cell": sub.index[mask],
        "cluster": final
    })

    out.to_csv(f"{name}_hdbscan_ab.tsv", sep="\t", index=False)

    plot_umap_3d(
        X2,
        final,
        f"{name} HDBSCAN + UMAP (scAB)",
        f"{name}_hdbscan_umap_ab.svg"
    )


# =========================
# RUN
# =========================
run_hdbscan("IVF-2-", "IVF-2C")
run_hdbscan("IVF-4-", "IVF-4C")