# python3

# ライブラリのインポート
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import tukey_hsd, mannwhitneyu
from statsmodels.stats.multicomp import pairwise_tukeyhsd  # 使えない
import japanize_matplotlib
import warnings
import pprint

warnings.simplefilter('ignore')
import sys

# snsのバージョンが0.13.0以降でないと、sns.boxplotのlinecolorの変数が指定不可
# print(sns.__version__)
# sys.exit()

"""パスたち"""
excel_path_mydata = "/Users/kotaro/PycharmProjects/iden/2023-11-02_実習B_1-6班.xls"
excel_path_reference = "/Users/kotaro/PycharmProjects/iden/実習B_qPCR_B3用参考データ.xls"
output_path = "/Users/kotaro/Desktop/遺伝/"
order = ["WT", "WT-PBS", "WT-Ecoli", "C", "C-PBS", "C-Ecoli"]

"""使う関数"""


def yuisa(x, threshold=[0.01, 0.05]):
    if x <= threshold[0]:
        return "**"
    elif x <= threshold[1]:
        return "*"
    else:
        return "n.s"


"""自班データの解析"""
print("【自班のデータ】")
# エクセルファイルを取得
df_mydata = pd.read_excel(excel_path_mydata, sheet_name="Results", skiprows=46, skipfooter=5)

# グループ情報の取得
df_mydata["char_Well Position"] = df_mydata["Well Position"].str[0]
lst_chars = list(df_mydata["char_Well Position"].unique())
lst_chars.sort()  # 念の為、ソート
lst_chars = lst_chars[:-1]  # H列の削除

dfs_mydata = []

# Relative Quantityの算出
for Group_char in lst_chars:
    df_extracted = df_mydata[df_mydata["Well Position"].str.contains(Group_char)][
        ["Sample Name", "Target Name", "Quantity"]]
    df_Dpt = df_extracted[df_extracted["Target Name"] == "Dpt"].rename(
        columns={"Target Name": "Target Name Dpt", "Quantity": "Quantity Dpt"})  # 各文字のDptの要素のdf
    df_pol = df_extracted[df_extracted["Target Name"] == "pol2"].rename(
        columns={"Target Name": "Target Name pol2", "Quantity": "Quantity pol2"})  # 各文字のpol2の要素のdf
    df_merged = pd.merge(df_Dpt, df_pol, on="Sample Name")  # "Sample Name"の下でマージ
    df_merged["Relative Quantity"] = df_merged["Quantity Dpt"] / df_merged["Quantity pol2"]
    dfs_mydata.append(df_merged)

# 結合したものをAから順にまとめたエクセルファイルの出力
df_mydata_concat = (pd.concat(dfs_mydata, axis=0))  # それぞれのdfをconcatenate
# df_mydata_concat.to_csv(output_path+"Relative_quantity.csv",encoding="cp932")

# サンプルごとにまとめたエクセルファイルの出力
df_mydata_samples = []
for i in df_mydata_concat["Sample Name"].unique():
    df_mydata_samples.append(df_mydata_concat[df_mydata_concat["Sample Name"] == i])

df_reference_concat_samples = (pd.concat(df_mydata_samples, axis=0))
df_reference_concat_samples.to_csv(output_path + "Relative_quantity(Sample_name).csv", encoding="cp932")

# 統計処理
group_mydata_sample_name = df_mydata_concat.groupby("Sample Name")

lst_mydata_stat = []
lst_mydata_stat_names = []

for i, v in group_mydata_sample_name:
    lst_mydata_stat_names.append(i)
    temp = group_mydata_sample_name.get_group(i)
    lst_mydata_stat.append(list(temp["Relative Quantity"]))

res = tukey_hsd(lst_mydata_stat[0], lst_mydata_stat[1], lst_mydata_stat[2], lst_mydata_stat[3], lst_mydata_stat[4],
                lst_mydata_stat[5])

# 統計結果の表示
df_mydata_pvalue = pd.DataFrame(data=res.pvalue, index=lst_mydata_stat_names, columns=lst_mydata_stat_names).reindex(
    index=order, columns=order)
for i in df_mydata_pvalue.columns: df_mydata_pvalue[i] = df_mydata_pvalue[i].apply(round, ndigits=4)
for i in df_mydata_pvalue.columns: df_mydata_pvalue[i] = df_mydata_pvalue[i].apply(yuisa)
print(df_mydata_pvalue)

# グラフの作成
# アスタリスクを入れる座標の取得
aster1_my = []
aster2_my = []
for i in range(len(df_mydata_pvalue)):
    for j in range(len(df_mydata_pvalue)):
        if df_mydata_pvalue.iat[i, j] == "*":
            if i > j:
                aster1_my.append([i, j])
        if df_mydata_pvalue.iat[i, j] == "**":
            if i > j:
                aster2_my.append([i, j])

sns.boxplot(data=df_mydata_concat, x="Sample Name", y="Relative Quantity",
            order=order, color='white', linecolor="black")
sns.swarmplot(data=df_mydata_concat, x="Sample Name", y="Relative Quantity",
              order=order, palette='Set2')

# 有意差の書き込み
cnt = 0.05
max_my = 27000
for key in aster1_my:
    plt.plot([key[0], key[1]], [max_my * (1 + cnt), max_my * (1 + cnt)], marker="|", color="black")
    plt.text((key[0] + key[1]) / 2, max_my * (1 + cnt - 0.025), "*", horizontalalignment="center",
             verticalalignment="bottom", fontsize="15")
    cnt += 0.05
for key in aster2_my:
    plt.plot([key[0], key[1]], [max_my * (1 + cnt), max_my * (1 + cnt)], marker="|", color="black")
    plt.text((key[0] + key[1]) / 2, max_my * (1 + cnt - 0.025), "**", horizontalalignment="center",
             verticalalignment="bottom", fontsize="15")
    cnt += 0.05

plt.title("各サンプルのmRNA量")
plt.text(4, -4000, "*:p<0.05, **:p<0.01", verticalalignment="top")
plt.tight_layout()
plt.savefig(output_path + "RT-RNA.jpg", dpi=300)
plt.show()
plt.close()

"""参考データの解析"""
print("【参考データ】")

# エクセルファイルを取得
df_reference = pd.read_excel(excel_path_reference, sheet_name="Results", skiprows=46, skipfooter=5)

# 使用していないグループの除去
df_reference = df_reference[~df_reference["Sample Name"].isna()]
df_reference = df_reference[["Well Position", "Sample Name", "Target Name", "Quantity"]]

# サンプル名称の統一
df_reference.loc[df_reference["Sample Name"] == "WT infection", "Sample Name"] = "WT-Ecoli"
df_reference.loc[df_reference["Sample Name"] == "WT PBS", "Sample Name"] = "WT-PBS"
df_reference.loc[df_reference["Sample Name"] == "C infection", "Sample Name"] = "C-Ecoli"
df_reference.loc[df_reference["Sample Name"] == "C PBS", "Sample Name"] = "C-PBS"

# Relative Quantityの導出
dfs_reference = []
for Group_char in list(df_reference["Sample Name"].unique()):
    df_extracted2 = df_reference[df_reference["Sample Name"] == Group_char][["Sample Name", "Target Name", "Quantity"]]
    df_Dpt = df_extracted2[df_extracted2["Target Name"] == "Dipt"].rename(
        columns={"Target Name": "Target Name Dpt"}).reset_index(drop=True)
    df_pol = df_extracted2[df_extracted2["Target Name"] == "pol2"].rename(
        columns={"Target Name": "Target Name pol2"}).reset_index(drop=True)
    for i in range(len(df_Dpt)):
        df_Dpt.loc[i, "Relative Quantity"] = df_Dpt.loc[i, "Quantity"] / df_pol.loc[i, "Quantity"]
    dfs_reference.append(df_Dpt)

df_reference_concat = pd.concat(dfs_reference, axis=0)

# 統計処理
lst_reference_stat = []
lst_reference_stat_names = []
group_reference_sample_name = df_reference_concat.groupby("Sample Name")
for i, v in group_reference_sample_name:
    temp_temp = []
    lst_reference_stat_names.append(i)
    for j in range(len(group_reference_sample_name.get_group(i))):
        temp_temp.append(group_reference_sample_name.get_group(i).loc[j, "Relative Quantity"])
    lst_reference_stat.append(temp_temp)
res = tukey_hsd(lst_reference_stat[0], lst_reference_stat[1], lst_reference_stat[2], lst_reference_stat[3],
                lst_reference_stat[4],
                lst_reference_stat[5])

# 統計処理の出力
df_reference_pvalue = pd.DataFrame(data=res.pvalue, index=lst_reference_stat_names,
                                   columns=lst_reference_stat_names).reindex(index=order, columns=order)
for i in df_reference_pvalue.columns: df_reference_pvalue[i] = df_reference_pvalue[i].apply(round, ndigits=4)
for i in df_reference_pvalue.columns: df_reference_pvalue[i] = df_reference_pvalue[i].apply(yuisa)
print(df_reference_pvalue)

# グラフの出力
aster1_ref = []
aster2_ref = []
for i in range(len(df_reference_pvalue)):
    for j in range(len(df_reference_pvalue)):
        if df_reference_pvalue.iat[i, j] == "*":
            if i > j:
                aster1_ref.append([i, j])
        if df_reference_pvalue.iat[i, j] == "**":
            if i > j:
                aster2_ref.append([i, j])

sns.boxplot(data=df_reference_concat, x="Sample Name", y="Relative Quantity",
            order=order, color='white', linecolor="black")
sns.swarmplot(data=df_reference_concat, x="Sample Name", y="Relative Quantity",
              order=order, palette='Set2')
cnt = 0.05
max_ref = 1.59
for key in aster1_ref:
    plt.plot([key[0], key[1]], [max_ref * (1 + cnt), max_ref * (1 + cnt)], marker="|", color="black")
    plt.text((key[0] + key[1]) / 2, max_ref * (1 + cnt - 0.025), "*", horizontalalignment="center",
             verticalalignment="bottom", fontsize="15")
    cnt += 0.05
for key in aster2_ref:
    plt.plot([key[0], key[1]], [max_ref * (1 + cnt), max_ref * (1 + cnt)], marker="|", color="black")
    plt.text((key[0] + key[1]) / 2, max_ref * (1 + cnt - 0.025), "**", horizontalalignment="center",
             verticalalignment="bottom", fontsize="15")
    cnt += 0.05

plt.title("各サンプルのmRNA量(参考データ)")
plt.yticks([0.2 * i for i in range(int(max_ref / 0.2) + 2)])
plt.text(4, -0.25, "*:p<0.05, **:p<0.01", verticalalignment="top")

plt.tight_layout()
plt.savefig(output_path + "RT-RNA(reference).jpg", dpi=300)
plt.show()
plt.close()
