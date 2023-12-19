import matplotlib.pyplot as plt
import japanize_matplotlib
import pandas as pd

output_path = "/Users/kotaro/PycharmProjects/iden/outputs/tables/"

LacZ_values = [37.523, 25.775, 27.368, 27.787, 18.733, 31.046]
LacZ_stages = [5, 8, 9, 11, 17, 13]
CG4319_values = [13.359, 18.704, 22.931]
CG4319_stages = [8, 11, 13]
CG12284_values = [19.545, 46.904, 37.888]
CG12284_stages = [8, 11, 5]

datas = {"LacZ": {"values": [37.523, 25.775, 27.368, 27.787, 18.733, 31.046], "stages": [5, 8, 9, 11, 17, 13]},
         "CG4319": {"values": [13.359, 18.704, 22.931], "stages": [8, 11, 13]},
         "CG12284": {"values": [19.545, 46.904, 37.888], "stages": [8, 11, 5]}
         }

def compare_same_stage():
    common_stages = set(LacZ_stages) & set(CG4319_stages) & set(CG12284_stages)
    for common_stage in common_stages:
        lacZ_index = LacZ_stages.index(common_stage)
        CG4319_index = CG4319_stages.index(common_stage)
        CG12284_index = CG12284_stages.index(common_stage)
        values = [LacZ_values[lacZ_index], CG4319_values[CG4319_index], CG12284_values[CG12284_index]]
        samples = ["LacZ", "CG4319", "CG12284"]
        plt.bar(samples, values)
        plt.title("各系統の顕微鏡画像の輝度(STAGE"+str(common_stage)+")")
        plt.xlabel("系統")
        plt.ylabel("輝度(a.u.)")
        plt.savefig(output_path + "graph_stage" + str(common_stage) + ".jpg", dpi=300)
        plt.show()
    return True



def print_table():
    df = pd.DataFrame(data=[[0 for i in range(6)] for j in range(3)], index=["LacZ", "CG4319", "CG12284"],
                      columns=["STAGE5", "STAGE8", "STAGE9", "STAGE11", "STAGE13", "STAGE17"])
    temp = {}
    for key in datas:
        temp[key] = []
        for i in sorted(datas[key]["stages"]):
            temp[key].append("STAGE"+str(i))

    for key in datas:
        values = datas[key]["values"]
        stages = datas[key]["stages"]
        for i in range(len(values)):
            value = values[i]
            stage = stages[i]
            df.loc[key, "STAGE"+str(stage)] = value

    for i in range(3):
        for j in range(6):
            if df.iloc[i,j] == 0:
                df.iloc[i, j] = "-"

    print(df,"\n")
    df_extract = df[list(set(temp["LacZ"])&set(temp["CG4319"])&set(temp["CG12284"]))]
    return df,df_extract


def Table(df,w=7,h=2,outputPath=output_path,title=""):
    fig, ax = plt.subplots(figsize=(w,h))
    # fig, ax = plt.subplots()
    ax.axis('off')
    ax.table(
        df.values,
        rowLabels=df.index,
        colLabels=df.columns,
        loc='center',
        bbox=[0.05,0.05,1,1]
    )
    plt.savefig(outputPath+"表"+title+".jpg", dpi=300)
    plt.show()


def graph():
    compare_same_stage()
    print_table()
    Table(print_table()[0])
    Table(print_table()[1], w=6, title="共通")

if __name__ == "__main__":
    graph()
