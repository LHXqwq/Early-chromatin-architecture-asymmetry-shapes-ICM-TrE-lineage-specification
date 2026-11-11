from collections import defaultdict, Counter
import matplotlib.pyplot as plt
from IPython.display import SVG, display

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['svg.hashsalt'] = 'hello'

input_file_name = "IVF-2C_hdbscan_ab.tsv"
output_file_name = "IVF_2C_AB_clustering_embryo_type_combinations_complete_only.txt"


# =========================
# 1. 读取并解析数据
# =========================
def read_and_parse(input_file_name):

    embryo_data = defaultdict(list)

    with open(input_file_name, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')

            # 跳过header
            if parts[0] == "ID":
                continue

            if len(parts) < 2:
                continue

            embryo_id = '-'.join(parts[0].split('-')[:3])

            # =========================
            # ✔ 关键修改：解析 ps1 / ps2
            # =========================
            cell_type = parts[1].split('-')[-1]

            embryo_data[embryo_id].append(cell_type)

    return embryo_data


# =========================
# 2. 只保留完整胚胎
# =========================
def count_complete_embryos(embryo_data):

    embryo_type_count = Counter()
    filtered_embryos = {}

    for embryo_id, cell_types in embryo_data.items():

        # 必须恰好两个细胞
        if len(cell_types) != 2:
            continue

        # 剔除 Undefined
        if "Undefined" in cell_types:
            continue

        filtered_embryos[embryo_id] = cell_types

        sorted_types = sorted(cell_types)

        if sorted_types[0] == sorted_types[1]:
            combination = f"{sorted_types[0]}-{sorted_types[0]}"
        else:
            combination = '-'.join(sorted_types)

        embryo_type_count[combination] += 1

    total_embryos = len(filtered_embryos)

    return filtered_embryos, embryo_type_count, total_embryos


# =========================
# 3. 写入结果文件
# =========================
def write_results(embryo_data, combination_counts, total_embryos, output_file_name):

    with open(output_file_name, 'w') as output_file:

        output_file.write("Complete embryos only (2 blastomeres, both passed QC)\n\n")

        for embryo_id, cell_types in embryo_data.items():

            sorted_types = sorted(cell_types)

            if sorted_types[0] == sorted_types[1]:
                combination = f"{sorted_types[0]}-{sorted_types[0]}"
            else:
                combination = '-'.join(sorted_types)

            output_file.write(f"{embryo_id}\t{combination}\n")

        output_file.write("\nType Combinations and Percentages:\n")

        for combination, count in combination_counts.items():
            percentage = (count / total_embryos) * 100
            output_file.write(f"{combination}: {count} ({percentage:.2f}%)\n")

    print("Results saved to:", output_file_name)


# =========================
# 4. 绘图
# =========================
def plot_pie_chart(combination_counts, total_embryos):

    labels = list(combination_counts.keys())
    sizes = list(combination_counts.values())

    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in range(len(labels))]

    plt.figure(figsize=(10, 8))

    explode = [0.05] * len(labels)

    wedges, texts, autotexts = plt.pie(
        sizes,
        labels=labels,
        autopct='%1.1f%%',
        startangle=140,
        colors=colors,
        explode=explode,
        textprops={'fontsize': 20}
    )

    # =========================
    # ✔ 只修改文字颜色
    # =========================
    for i, text in enumerate(texts):
        if labels[i] == "AB1-AB2":
            text.set_color('red')
            text.set_fontweight('bold')

    plt.axis('equal')

    plt.title(
        f'IVF 2C Complete Embryos (n=29)\nCell Type Combinations (scAB) hierarchical',
        fontsize=16,
        pad=20
    )
    plt.savefig("ivf2c_scAB_hdbscan.svg", dpi=300, bbox_inches='tight', transparent=True)
    plt.show()


# =========================
# 5. 主函数
# =========================
def main():

    embryo_data = read_and_parse(input_file_name)

    filtered_embryos, combination_counts, total_embryos = \
        count_complete_embryos(embryo_data)

    print("Number of complete embryos:", total_embryos)

    write_results(
        filtered_embryos,
        combination_counts,
        total_embryos,
        output_file_name
    )

    plot_pie_chart(combination_counts, total_embryos)


if __name__ == "__main__":
    main()

from collections import defaultdict, Counter
import matplotlib.pyplot as plt

input_file_name = "IVF-4C_hdbscan_ab.tsv"
output_file_name = "IVF_4C_AB_clustering_counts_complete_only.txt"


# =========================
# 1. 读取数据
# =========================
def read_and_parse(input_file_name):

    embryo_data = defaultdict(list)

    with open(input_file_name, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')

            if parts[0] == "ID":
                continue

            if len(parts) < 2:
                continue

            embryo_id = '-'.join(parts[0].split('-')[:3])

            # ✔ ps1 / ps2 / Undefined
            cell_type = parts[1].split('-')[-1]

            embryo_data[embryo_id].append(cell_type)

    return embryo_data


# =========================
# 2. 只保留完整 4C
# =========================
def filter_embryos(embryo_data):

    filtered = {}
    counts = Counter()

    for eid, types in embryo_data.items():

        if len(types) != 4:
            continue

        if "Undefined" in types:
            continue

        filtered[eid] = types

        c = Counter(types)

        ps1 = c.get("AB1", 0)
        ps2 = c.get("AB2", 0)

        label = f"AB1({ps1})-AB2({ps2})"

        counts[label] += 1

    return filtered, counts, len(filtered)


# =========================
# 3. 输出文件
# =========================
def write_results(filtered, counts, total, out_file):

    with open(out_file, 'w') as f:

        f.write("IVF 4C embryos (ps counts only)\n\n")

        for eid, types in filtered.items():

            c = Counter(types)
            ps1 = c.get("AB1", 0)
            ps2 = c.get("AB2", 0)

            f.write(f"{eid}\tAB1({ps1})-AB2({ps2})\n")

        f.write("\nSummary:\n")

        for k, v in counts.items():
            f.write(f"{k}: {v} ({v/total*100:.2f}%)\n")

    print("Saved:", out_file)


# =========================
# 4. 饼图
# =========================
def plot_pie(counts, total):

    labels = list(counts.keys())
    sizes = list(counts.values())

    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in range(len(labels))]

    plt.figure(figsize=(10, 8))

    explode = [0.03] * len(labels)

    wedges, texts, autotexts = plt.pie(
        sizes,
        labels=labels,
        autopct='%1.1f%%',
        startangle=140,
        colors=colors,
        explode=explode,
        textprops={'fontsize': 20}
    )


    # =========================
    # ✔ 只要同时存在 ps1 + ps2 就标红
    # =========================
    for i, text in enumerate(texts):

        label = labels[i]

        try:
            part1, part2 = label.split('-')

            ps1 = int(part1.split('(')[1].split(')')[0])
            ps2 = int(part2.split('(')[1].split(')')[0])

            # ✔ 关键逻辑
            if ps1 > 0 and ps2 > 0:
                text.set_color('red')
                text.set_fontweight('bold')

        except:
            continue

    plt.axis('equal')

    plt.title(f"IVF 4C Complete Embryos (n=14)\nCell Type Combinations (scAB) hierarchical", fontsize=16,pad=20)
    plt.savefig("ivf4c_scAB_hdbscan.svg", dpi=300, bbox_inches='tight', transparent=True)
    plt.show()


# =========================
# 5. main
# =========================
def main():

    data = read_and_parse(input_file_name)

    filtered, counts, total = filter_embryos(data)

    print("Complete embryos:", total)

    write_results(filtered, counts, total, output_file_name)

    plot_pie(counts, total)


if __name__ == "__main__":
    main()