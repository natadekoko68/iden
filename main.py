# python3

# ライブラリのインポート
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import tukey_hsd
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import japanize_matplotlib
import warnings
warnings.simplefilter('ignore')
import sys

# snsのバージョンが0.13.0以降でないと、sns.boxplotのlinecolorの変数が指定不可
# print(sns.__version__)
# sys.exit()

"""自班データの解析"""
# エクセルファイルを取得
excel_path_jihan = "/Users/kotaro/PycharmProjects/iden/2023-11-02_実習B_1-6班.xls"
df_jihan = pd.read_excel(excel_path_jihan, sheet_name="Results", skiprows=46, skipfooter=5)

# グループ情報の取得
df_jihan["char_Well Position"] = df_jihan["Well Position"].str[0]
lst_chars = list(df_jihan["char_Well Position"].unique())
lst_chars.sort()  # 念の為、ソート
lst_chars = lst_chars[:-1]  # H列の削除

dfs_jihan = []

# Relative Quantityの算出
for Group_char in lst_chars:
    df_extracted = df_jihan[df_jihan["Well Position"].str.contains(Group_char)][["Sample Name", "Target Name", "Quantity"]]
    df_Dpt = df_extracted[df_extracted["Target Name"] == "Dpt"].rename(
        columns={"Target Name": "Target Name Dpt", "Quantity": "Quantity Dpt"})  # 各文字のDptの要素のdf
    df_pol = df_extracted[df_extracted["Target Name"] == "pol2"].rename(
        columns={"Target Name": "Target Name pol2", "Quantity": "Quantity pol2"})  # 各文字のpol2の要素のdf
    df_merged = pd.merge(df_Dpt, df_pol, on="Sample Name")  # "Sample Name"の下でマージ
    df_merged["Relative Quantity"] = df_merged["Quantity Dpt"] / df_merged["Quantity pol2"]
    dfs_jihan.append(df_merged)

# 結合したものをAから順にまとめたエクセルファイルの出力
df_jihan_concated = (pd.concat(dfs_jihan, axis=0))  # それぞれのdfをconcatenate
# df_jihan_concated.to_csv("/Users/kotaro/Desktop/遺伝R/Relative_quantity.csv",encoding="cp932")

# サンプルごとにまとめたエクセルファイルの出力
df_jihan_samples = []
for i in df_jihan_concated["Sample Name"].unique():
    df_jihan_samples.append(df_jihan_concated[df_jihan_concated["Sample Name"] == i])

df_sanko_concat_samples = (pd.concat(df_jihan_samples, axis=0))
df_sanko_concat_samples.to_csv("/Users/kotaro/Desktop/遺伝/Relative_quantity(Sample_name).csv", encoding="cp932")

# グラフの作成
sns.boxplot(data=df_jihan_concated, x="Sample Name", y="Relative Quantity",
            order=["WT", "WT-Ecoli", "WT-PBS", "C", "C-Ecoli", "C-PBS"], color='white', linecolor="black")
sns.swarmplot(data=df_jihan_concated, x="Sample Name", y="Relative Quantity",
              order=["WT", "WT-Ecoli", "WT-PBS", "C", "C-Ecoli", "C-PBS"], palette='Set2')
plt.title("各サンプルのmRNA量")
for i in range(1, 6):
    plt.text(i,
             df_jihan_concated[df_jihan_concated["Sample Name"] == list(["WT", "WT-Ecoli", "WT-PBS", "C", "C-Ecoli", "C-PBS"])[i]][
                 "Relative Quantity"].max() * 1.04, " n.s.", verticalalignment='bottom', horizontalalignment="center")
plt.ylim([-1000, 30000])
plt.tight_layout()

plt.savefig("/Users/kotaro/Desktop/遺伝/RT-RNA Result.jpg")
plt.show()
plt.close()

# 統計処理
group_jihan_sample_name = df_jihan_concated.groupby("Sample Name")

lst_jihan_stat = []
lst_jihan_stat_names = []

for i, v in group_jihan_sample_name:
    # print(i)
    lst_jihan_stat_names.append(i)
    temp = group_jihan_sample_name.get_group(i)
    lst_jihan_stat.append(list(temp["Relative Quantity"]))

res = tukey_hsd(lst_jihan_stat[0], lst_jihan_stat[1], lst_jihan_stat[2], lst_jihan_stat[3], lst_jihan_stat[4], lst_jihan_stat[5])
print("自班のデータ")
print("番号の対応 : ", lst_jihan_stat_names)
# res.confidence_interval(confidence_level=0.99)
print(res)

"""参考データの解析"""
# エクセルファイルを取得
excel_path_sanko = "/Users/kotaro/PycharmProjects/iden/実習B_qPCR_B3用参考データ.xls"
df_sanko = pd.read_excel(excel_path_sanko, sheet_name="Results", skiprows=46, skipfooter=5)

# 使用していないグループの除去
df_sanko = df_sanko[~df_sanko["Sample Name"].isna()]
df_sanko = df_sanko[["Well Position", "Sample Name", "Target Name", "Quantity"]]

# サンプル名称の統一
df_sanko.loc[df_sanko["Sample Name"] == "WT infection", "Sample Name"] = "WT-Ecoli"
df_sanko.loc[df_sanko["Sample Name"] == "WT PBS", "Sample Name"] = "WT-PBS"
df_sanko.loc[df_sanko["Sample Name"] == "C infection", "Sample Name"] = "C-Ecoli"
df_sanko.loc[df_sanko["Sample Name"] == "C PBS", "Sample Name"] = "C-PBS"

# Relative Quantityの導出
dfs_sanko = []
for Group_char in list(df_sanko["Sample Name"].unique()):
    df_extracted2 = df_sanko[df_sanko["Sample Name"] == Group_char][["Sample Name", "Target Name", "Quantity"]]
    df_Dpt = df_extracted2[df_extracted2["Target Name"] == "Dipt"].rename(
        columns={"Target Name": "Target Name Dpt"}).reset_index(drop=True)
    df_pol = df_extracted2[df_extracted2["Target Name"] == "pol2"].rename(
        columns={"Target Name": "Target Name pol2"}).reset_index(drop=True)
    for i in range(len(df_Dpt)):
        df_Dpt.loc[i, "Relative Quantity"] = df_Dpt.loc[i, "Quantity"] / df_pol.loc[i, "Quantity"]
    dfs_sanko.append(df_Dpt)

df_sanko_concat = pd.concat(dfs_sanko, axis=0)

# グラフの描画
sns.boxplot(data=df_sanko_concat, x="Sample Name", y="Relative Quantity",
            order=["WT", "WT-PBS", "WT-Ecoli", "C", "C-PBS", "C-Ecoli"], color='white', linecolor="black")
sns.swarmplot(data=df_sanko_concat, x="Sample Name", y="Relative Quantity",
              order=["WT", "WT-PBS", "WT-Ecoli", "C", "C-PBS", "C-Ecoli"], palette='Set2')

plt.title("各サンプルのmRNA量(参考データ)")
plt.tight_layout()
plt.savefig("/Users/kotaro/Desktop/遺伝/RT-RNA(reference).jpg")
plt.show()
plt.close()

# 統計処理
lst_sanko_stat = []
lst_sanko_stat_names = []
group_sanko_sample_name = df_sanko_concat.groupby("Sample Name")
for i, v in group_sanko_sample_name:
    temp_temp = []
    lst_sanko_stat_names.append(i)
    for j in range(len(group_sanko_sample_name.get_group(i))):
        temp_temp.append(group_sanko_sample_name.get_group(i).loc[j, "Relative Quantity"])
    lst_sanko_stat.append(temp_temp)
res = tukey_hsd(lst_sanko_stat[0], lst_sanko_stat[1], lst_sanko_stat[2], lst_sanko_stat[3], lst_sanko_stat[4], lst_sanko_stat[5])
# res.confidence_interval(confidence_level=0.99)
print("参考データ")
print("番号の対応 : ", lst_sanko_stat_names)
print(res)
