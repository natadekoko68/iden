# python3

# ライブラリのインポート
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import tukey_hsd
import japanize_matplotlib
import warnings

warnings.simplefilter('ignore')

# version check: snsのバージョンが0.13.0以降でないと、sns.boxplotのlinecolorの変数が指定不可
assert ((int(sns.__version__.split(".")[1]) >= 13) or (int(sns.__version__.split(".")[0]) >= 1))

"""パスたち"""
excel_path_mydata = "/Users/kotaro/PycharmProjects/iden/inputs/2023-11-02_実習B_1-6班.xls"
excel_path_reference = "/Users/kotaro/PycharmProjects/iden/inputs/実習B_qPCR_B3用参考データ.xls"
output_path = "/Users/kotaro/PycharmProjects/iden/outputs/tables/"
order = ["WT", "WT-PBS", "WT-Ecoli", "C", "C-PBS", "C-Ecoli"]

"""使う関数"""


def yuisa(x, threshold=[0.01, 0.05]):
    if x <= threshold[0]:
        return "**"
    elif x <= threshold[1]:
        return "*"
    else:
        return "n.s"


def stat_display(df, order=order):
    """
    :param df: df もとデータ
    :param order: list df_pvalueで表示する順序
    :return:
    """
    group = df.groupby("Sample Name")
    lst_stat = []
    lst_stat_names = []
    for i, v in group:
        lst_stat_names.append(i)
        lst_stat.append(list(group.get_group(i)["Relative Quantity"]))
    res = tukey_hsd(lst_stat[0], lst_stat[1], lst_stat[2], lst_stat[3], lst_stat[4], lst_stat[5])
    df_pvalue = pd.DataFrame(data=res.pvalue, index=lst_stat_names, columns=lst_stat_names).reindex(index=order,
                                                                                                    columns=order)
    for i in df_pvalue.columns: df_pvalue[i] = df_pvalue[i].apply(round, ndigits=4)
    for i in df_pvalue.columns: df_pvalue[i] = df_pvalue[i].apply(yuisa)
    print(df_pvalue)
    aster1 = []
    aster2 = []
    for i in range(len(df_pvalue)):
        for j in range(len(df_pvalue)):
            if df_pvalue.iat[i, j] == "*":
                if i > j:
                    aster1.append([i, j])
            if df_pvalue.iat[i, j] == "**":
                if i > j:
                    aster2.append([i, j])
    return df_pvalue, aster1, aster2


def graph(data, aster1, aster2, cnt=0.05, df_max=27000, txt_pos_y=-4000, title="各サンプルのmRNA量",
          output_path=output_path, add_text="", order=order):
    """
    :param data: df　出力するグラフのもとデータ
    :param aster1: list *を描画する群のデータ
    :param aster2: list **を描画する群のデータ
    :param cnt: float 平行線をずらして描画する程度
    :param df_max: float dfの最大値
    :param txt_pos_y: float "*:p<0.05, **:p<0.01"を入れる位置のy座標
    :param title: str グラフに表示するタイトル
    :param output_path: str 出力ディレクトリ
    :param add_text: str 保存するときに付加する文章
    :param order: グラフに表示する順序
    """
    for aster1_pos in aster1:
        plt.plot([aster1_pos[0], aster1_pos[1]], [df_max * (1 + cnt), df_max * (1 + cnt)], marker="|", color="black")
        plt.text((aster1_pos[0] + aster1_pos[1]) / 2, df_max * (1 + cnt - 0.025), "*", horizontalalignment="center",
                 verticalalignment="bottom", fontsize="15")
        cnt += 0.05
    for aster2_pos in aster2:
        plt.plot([aster2_pos[0], aster2_pos[1]], [df_max * (1 + cnt), df_max * (1 + cnt)], marker="|", color="black")
        plt.text((aster2_pos[0] + aster2_pos[1]) / 2, df_max * (1 + cnt - 0.025), "**", horizontalalignment="center",
                 verticalalignment="bottom", fontsize="15")
        cnt += 0.05
    sns.boxplot(data=data, x="Sample Name", y="Relative Quantity",
                order=order, color='white', linecolor="black")
    sns.swarmplot(data=data, x="Sample Name", y="Relative Quantity",
                  order=order, palette='Set2')
    plt.title(title)
    plt.text(4, txt_pos_y, "*:p<0.05, **:p<0.01", verticalalignment="top")
    plt.tight_layout()
    plt.savefig(output_path + "RT-RNA" + str(add_text) + ".jpg", dpi=300)
    plt.show()
    plt.close()


def get_mydata():
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
    df_mydata_samples = []
    for i in df_mydata_concat["Sample Name"].unique():
        q = df_mydata_concat[df_mydata_concat["Sample Name"] == i]["Relative Quantity"].quantile(0.99)
        df_mydata_samples.append(df_mydata_concat[df_mydata_concat["Sample Name"] == i][
                                     df_mydata_concat[df_mydata_concat["Sample Name"] == i]["Relative Quantity"] <= q])
    return df_mydata_concat


def heikin_graph(data, outlier=False, quantile=0.99, add_text=""):
    df_mydata_concat = get_mydata()
    temp = []
    for i in order:
        if outlier:
            temp.append(df_mydata_concat[df_mydata_concat["Sample Name"] == i][
                            df_mydata_concat[df_mydata_concat["Sample Name"] == i]["Relative Quantity"] <=
                            df_mydata_concat[df_mydata_concat["Sample Name"] == i]["Relative Quantity"].quantile(
                                quantile)]["Relative Quantity"].mean())
        else:
            temp.append(data.groupby("Sample Name").get_group(i)["Relative Quantity"].mean())
    plt.bar(order, temp)
    plt.title("各サンプルの相対的RNA発現量の平均値" + add_text)
    plt.xlabel("系統")
    plt.ylabel("相対的RNA発現量(a.u.)")
    if outlier:
        plt.savefig(output_path + "平均値(外れ値除外)" + add_text + ".jpg", dpi=300)
    else:
        plt.savefig(output_path + "平均値" + add_text + ".jpg", dpi=300)
    plt.show()
    return True


def get_referencedata():
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
        df_extracted2 = df_reference[df_reference["Sample Name"] == Group_char][
            ["Sample Name", "Target Name", "Quantity"]]
        df_Dpt = df_extracted2[df_extracted2["Target Name"] == "Dipt"].rename(
            columns={"Target Name": "Target Name Dpt"}).reset_index(drop=True)
        df_pol = df_extracted2[df_extracted2["Target Name"] == "pol2"].rename(
            columns={"Target Name": "Target Name pol2"}).reset_index(drop=True)
        for i in range(len(df_Dpt)):
            df_Dpt.loc[i, "Relative Quantity"] = df_Dpt.loc[i, "Quantity"] / df_pol.loc[i, "Quantity"]
        dfs_reference.append(df_Dpt)
        df_reference_concat = pd.concat(dfs_reference, axis=0)
    return df_reference_concat

def processing():
    df_mydata_concat = get_mydata()
    stat_display(df_mydata_concat)
    print("")
    heikin_graph(df_mydata_concat, outlier=True)
    print("")
    heikin_graph(df_mydata_concat)
    print("")

    df_reference_concat = get_referencedata()
    stat_display(df_reference_concat)
    print("")
    heikin_graph(df_reference_concat, outlier=True, add_text="(参考データ)")
    print("")
    heikin_graph(df_reference_concat, add_text="(参考データ)")
    print("")


if __name__ == "__main__":
    processing()