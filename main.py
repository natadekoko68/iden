# ライブラリのインポート
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import tukey_hsd
import japanize_matplotlib


"""自班データの解析"""
# エクセルファイルを取得
excel_path = "/Users/kotaro/Desktop/遺伝R/2023-11-02_実習B_1-6班.xls"
df3 = pd.read_excel(excel_path,sheet_name="Results",skiprows=46,skipfooter=5)

# グループ情報の取得
df3["char_Well Position"] = df3["Well Position"].str[0]
chars = list(df3["char_Well Position"].unique())
chars.sort() #念の為、ソート
chars = chars[:-1] #H列の削除

dfs = []

# Relative Quantityの算出
for Group_char in chars:
    df_extracted = df3[df3["Well Position"].str.contains(Group_char)][["Sample Name", "Target Name", "Quantity"]]
    df_Dpt = df_extracted[df_extracted["Target Name"] == "Dpt"].rename(
        columns={"Target Name": "Target Name Dpt", "Quantity": "Quantity Dpt"}) #各文字のDptの要素のdf
    df_pol = df_extracted[df_extracted["Target Name"] == "pol2"].rename(
        columns={"Target Name": "Target Name pol2", "Quantity": "Quantity pol2"}) #各文字のpol2の要素のdf
    df_merged = pd.merge(df_Dpt, df_pol, on="Sample Name") #"Sample Name"の下でマージ
    df_merged["Relative Quantity"] = df_merged["Quantity Dpt"] / df_merged["Quantity pol2"]
    dfs.append(df_merged)

# 結合したものをAから順にまとめたエクセルファイルの出力
df_concated = (pd.concat(dfs, axis=0)) #それぞれのdfをconcatenate
# df_concated.to_csv("/Users/kotaro/Desktop/遺伝R/Relative_quantity.csv",encoding="cp932")

# サンプルごとにまとめたエクセルファイルの出力
dfs3 = []
for i in df_concated["Sample Name"].unique():
    dfs3.append(df_concated[df_concated["Sample Name"] == i])

df_concated3 = (pd.concat(dfs3, axis=0))
df_concated3.to_csv("/Users/kotaro/Desktop/遺伝R/Relative_quantity(Sample_name).csv",encoding="cp932")


# グラフの作成
sns.boxplot(data=df_concated, x="Sample Name", y="Relative Quantity", order=["WT","WT-Ecoli","WT-PBS","C","C-Ecoli","C-PBS"],color='white')
sns.swarmplot(data=df_concated, x="Sample Name", y="Relative Quantity", order=["WT","WT-Ecoli","WT-PBS","C","C-Ecoli","C-PBS"], palette='Set2')
plt.title("各サンプルのmRNA量")
for i in range(1,6):
    plt.text(i,df_concated[df_concated["Sample Name"] == list(["WT","WT-Ecoli","WT-PBS","C","C-Ecoli","C-PBS"])[i]]["Relative Quantity"].max()*1.04," n.s.",verticalalignment='bottom',horizontalalignment="center")
plt.ylim([-1000,30000])
plt.tight_layout()

plt.savefig("/Users/kotaro/Desktop/遺伝R/RT-RNA Result.jpg")
plt.show()
plt.close()


# 統計処理
group = df_concated.groupby("Sample Name")

temps = []
temps_name = []

for i,v in group:
    # print(i)
    temps_name.append(i)
    temp = group.get_group(i)
    temps.append(list(temp["Relative Quantity"]))


res = tukey_hsd(temps[0],temps[1],temps[2],temps[3],temps[4],temps[5])
print("番号の対応 : ",temps_name)
# print(res.confidence_interval(confidence_level=0.99))
print(res)

"""参考データの解析"""
# エクセルファイルを取得
excel_path = "/Users/kotaro/Downloads/実習B_qPCR_B3用参考データ (1).xls"
df_temp = pd.read_excel(excel_path, sheet_name="Results", skiprows=46, skipfooter=5)

# 使用していないグループの除去
df_temp = df_temp[~df_temp["Sample Name"].isna()]
df_temp = df_temp[["Well Position","Sample Name", "Target Name", "Quantity"]]

# サンプル名称の統一
df_temp.loc[df_temp["Sample Name"] == "WT infection","Sample Name"] = "WT-Ecoli"
df_temp.loc[df_temp["Sample Name"] == "WT PBS","Sample Name"] = "WT-PBS"
df_temp.loc[df_temp["Sample Name"] == "C infection","Sample Name"] = "C-Ecoli"
df_temp.loc[df_temp["Sample Name"] == "C PBS","Sample Name"] = "C-PBS"

# Relative Quantityの導出
dfs2 = []
for Group_char in list(df_temp["Sample Name"].unique()):
    df_extracted2 = df_temp[df_temp["Sample Name"] == Group_char][["Sample Name", "Target Name", "Quantity"]]
    df_Dpt = df_extracted2[df_extracted2["Target Name"] == "Dipt"].rename(
        columns={"Target Name": "Target Name Dpt"}).reset_index(drop=True)
    df_pol = df_extracted2[df_extracted2["Target Name"] == "pol2"].rename(
        columns={"Target Name": "Target Name pol2"}).reset_index(drop=True)
    for i in range(len(df_Dpt)):
        df_Dpt.loc[i,"Relative Quantity"] = df_Dpt.loc[i,"Quantity"]/df_pol.loc[i,"Quantity"]
    dfs2.append(df_Dpt)

df_concated2 = pd.concat(dfs2, axis=0)

# グラフの描画
sns.boxplot(data=df_concated2, x="Sample Name", y="Relative Quantity", order=["WT","WT-PBS","WT-Ecoli","C","C-PBS","C-Ecoli"],color='white')
sns.swarmplot(data=df_concated2, x="Sample Name", y="Relative Quantity", order=["WT","WT-PBS","WT-Ecoli","C","C-PBS","C-Ecoli"], palette='Set2')
plt.tight_layout()
plt.title("各サンプルのmRNA量(参考データ)")
plt.savefig("/Users/kotaro/Desktop/遺伝R/RT-RNA(reference).jpg")
plt.show()
plt.close()

# 統計処理
temps2 = []
group_temp = df_concated2.groupby("Sample Name")
for i,v in group_temp:
    temp_temp = []
    print(i)
    for j in range(len(group_temp.get_group(i))):
        temp_temp.append(group_temp.get_group(i).loc[j,"Relative Quantity"])
    temps2.append(temp_temp)
res = tukey_hsd(temps2[0],temps2[1],temps2[2],temps2[3],temps2[4],temps2[5])
# print(res.confidence_interval(confidence_level=0.99))
print(res)


