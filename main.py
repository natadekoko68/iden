
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scipy.stats import tukey_hsd
import japanize_matplotlib

#     temps.append(list(temp["Value"]))

# res = tukey_hsd(temps[0],temps[1],temps[2],temps[3])
# print(res.confidence_interval(confidence_level=0.99))
# print(res)



excel_path = "/Users/kotaro/Desktop/遺伝R/2023-11-02_実習B_1-6班.xls"
# excel_path = "/Users/kotaro/Downloads/実習B_qPCR_B3用参考データ (1).xls"
# df1 = pd.read_excel(excel_path,sheet_name="Sample Setup",skiprows=46)
# df2 = pd.read_excel(excel_path,sheet_name="Amplification Data",skiprows=46)
df3 = pd.read_excel(excel_path,sheet_name="Results",skiprows=46,skipfooter=5)
# df4 = pd.read_excel(excel_path,sheet_name="Melt Curve Raw Data",skiprows=46)




# print(df1)
# print(df1.columns)

# print(df2)
# print(df2.columns)

df3["char_Well Position"] = df3["Well Position"].str[0]
# print(df3)
chars = list(df3["char_Well Position"].unique())
# chars = set()
# for i in range(len(df3)):
#     chars.add(str(df3.loc[i,"Well Position"])[0])
# chars = list(chars)
chars.sort()
chars = chars[:-1]

dfs = []

# print(chars)

for Group_char in chars:
    # print(Group_char)
    df_extracted = df3[df3["Well Position"].str.contains(Group_char)][["Sample Name", "Target Name", "Quantity"]]
    df_Dpt = df_extracted[df_extracted["Target Name"] == "Dpt"].rename(
        columns={"Target Name": "Target Name Dpt", "Quantity": "Quantity Dpt"})
    df_pol = df_extracted[df_extracted["Target Name"] == "pol2"].rename(
        columns={"Target Name": "Target Name pol2", "Quantity": "Quantity pol2"})
    # print(df_Dpt)
    # print(df_pol)

    df_merged = pd.merge(df_Dpt, df_pol, on="Sample Name")
    df_merged["Relative Quantity"] = df_merged["Quantity Dpt"] / df_merged["Quantity pol2"]
    # print(df_merged)
    dfs.append(df_merged)
    dfs

df_concated = (pd.concat(dfs, axis=0))

dfs3 = []
for i in df_concated["Sample Name"].unique():
    dfs3.append(df_concated[df_concated["Sample Name"]==i])

df_concated3 = (pd.concat(dfs3, axis=0))
df_concated3.to_csv("/Users/kotaro/Desktop/遺伝R/Relative quantity3.csv",encoding="cp932")



for key in df_concated["Sample Name"].unique():
    print(df_concated[df_concated["Sample Name"]==key])

df_concated.to_csv("/Users/kotaro/Desktop/遺伝R/Relative quantity.csv",encoding="cp932")


# sns.violinplot(data=df_concated, x="Sample Name", y="Relative Quantity")#, hue="alive", fill=False)
# plt.savefig("/Users/kotaro/Desktop/遺伝R/violin.jpg")
# plt.close()
# plt.show()

# sns.catplot(data=df_concated, x="Sample Name", y="Relative Quantity")#, hue="smoker")
# plt.savefig("/Users/kotaro/Desktop/遺伝R/swarm.jpg")
# plt.close()
# plt.show()


sns.boxplot(data=df_concated, x="Sample Name", y="Relative Quantity", order=["WT","WT-Ecoli","WT-PBS","C","C-Ecoli","C-PBS"],color='white')
sns.swarmplot(data=df_concated, x="Sample Name", y="Relative Quantity", order=["WT","WT-Ecoli","WT-PBS","C","C-Ecoli","C-PBS"], palette='Set2')
# plt.legend()
plt.title("各サンプルのmRNA量")
for i in range(1,6):
    plt.text(i,df_concated[df_concated["Sample Name"] == list(["WT","WT-Ecoli","WT-PBS","C","C-Ecoli","C-PBS"])[i]]["Relative Quantity"].max()*1.04," n.s.",verticalalignment='bottom',horizontalalignment="center")
plt.tight_layout()
plt.ylim([-1000,30000])
plt.savefig("/Users/kotaro/Desktop/遺伝R/Sample.jpg")
plt.show()
plt.close()

# fig = plt.figure()
# sns.boxplot(data=df_concated, x="Sample Name", y="Relative Quantity", inner="quartile", color='white')
# sns.swarmplot(data=df_concated, x="Sample Name", y="Relative Quantity", palette="Set2")
# plt.tight_layout()
# plt.show()

"""
new_dfs = []
samples = df_concated["Sample Name"].unique()
for sample in samples:
    temp = df_concated.groupby("Sample Name").get_group(sample)
    q_99 = temp["Relative Quantity"].quantile(0.999)
    q_01 = temp["Relative Quantity"].quantile(0.001)
    new_dfs.append(temp[(temp["Relative Quantity"] <= q_99)&(temp["Relative Quantity"] >= q_01)])
print(new_dfs)

df_concated_new = pd.concat(new_dfs, axis=0)
sns.violinplot(data=df_concated_new, x="Sample Name", y="Relative Quantity")#, hue="alive", fill=False)
# plt.savefig("/Users/kotaro/Desktop/遺伝R/violin_外れ値除外.jpg")
plt.close()
# plt.show()

sns.catplot(data=df_concated_new, x="Sample Name", y="Relative Quantity")
# plt.savefig("/Users/kotaro/Desktop/遺伝R/swarm_外れ値除外.jpg")
plt.close()
# plt.show()
"""


group = df_concated.groupby("Sample Name")

temps = []
temps_name = []

for i,v in group:
    # print(i)
    temps_name.append(i)
    temp = group.get_group(i)
    temps.append(list(temp["Relative Quantity"]))


res = tukey_hsd(temps[0],temps[1],temps[2],temps[3],temps[4],temps[5])
# print(res.confidence_interval(confidence_level=0.99))
print(res)
# res_df = pd.DataFrame(res)
# print(res_df)


"""参考データの解析"""

excel_path = "/Users/kotaro/Downloads/実習B_qPCR_B3用参考データ (1).xls"

df_temp = pd.read_excel(excel_path, sheet_name="Results", skiprows=46, skipfooter=5)
# print(df_temp[["Sample Name", "Target Name", "Quantity"]])
df_temp = df_temp[~df_temp["Sample Name"].isna()]
df_temp = df_temp[["Well Position","Sample Name", "Target Name", "Quantity"]]
# ["WT","WT-Ecoli","WT-PBS","C","C-Ecoli","C-PBS"]

df_temp.loc[df_temp["Sample Name"] == "WT infection","Sample Name"] = "WT-Ecoli"
df_temp.loc[df_temp["Sample Name"] == "WT PBS","Sample Name"] = "WT-PBS"
df_temp.loc[df_temp["Sample Name"] == "C infection","Sample Name"] = "C-Ecoli"
df_temp.loc[df_temp["Sample Name"] == "C PBS","Sample Name"] = "C-PBS"

#
# df_temp["char_Well Position"] = df_temp["Well Position"].str[0]
chars = list(df_temp["Sample Name"].unique())
chars.sort()
print(chars)

dfs2 = []
for Group_char in chars:
    df_extracted2 = df_temp[df_temp["Sample Name"] == Group_char][["Sample Name", "Target Name", "Quantity"]]
    # print(df_extracted2)
    df_Dpt = df_extracted2[df_extracted2["Target Name"] == "Dipt"].rename(
        columns={"Target Name": "Target Name Dpt"}).reset_index(drop=True)
    df_pol = df_extracted2[df_extracted2["Target Name"] == "pol2"].rename(
        columns={"Target Name": "Target Name pol2"}).reset_index(drop=True)
    # print(df_Dpt)
    for i in range(len(df_Dpt)):
        df_Dpt.loc[i,"Relative Quantity"] = df_Dpt.loc[i,"Quantity"]/df_pol.loc[i,"Quantity"]
    dfs2.append(df_Dpt)
# print(dfs2)

df_concated2 = pd.concat(dfs2, axis=0)

# for i in range(1,6):
#     plt.text(i,df_concated[df_concated["Sample Name"] == list(["WT","WT-Ecoli","WT-PBS","C","C-Ecoli","C-PBS"])[i]]["Relative Quantity"].max()*1.04," n.s.",verticalalignment='bottom',horizontalalignment="center")

sns.boxplot(data=df_concated2, x="Sample Name", y="Relative Quantity", order=["WT","WT-PBS","WT-Ecoli","C","C-PBS","C-Ecoli"],color='white') #, order=["WT","WT-PBS","WT-Ecoli","C","C-PBS","C-Ecoli"]
sns.swarmplot(data=df_concated2, x="Sample Name", y="Relative Quantity", order=["WT","WT-PBS","WT-Ecoli","C","C-PBS","C-Ecoli"], palette='Set2')
plt.tight_layout()
plt.show()
# plt.close()

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
print(temps2)

print(df_concated2)

df_concated2.to_csv("/Users/kotaro/Desktop/csv.csv")


