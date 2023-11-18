import matplotlib.pyplot as plt
import japanize_matplotlib

output_path = "/Users/kotaro/Desktop/"

LacZ_values = [37.523, 25.775, 27.368, 27.787, 18.733, 31.046]
LacZ_stages = [5, 8, 9, 11, 17, 13]
CG4319_values = [13.359, 18.704, 22.931]
CG4319_stages = [8, 11, 13]
CG12284_values = [19.545, 46.904, 37.888]
CG12284_stages = [8, 11, 5]

common_stages = list(set(LacZ_stages) & set(CG4319_stages) & set(CG12284_stages))
for common_stage in common_stages:
    lacZ_index = LacZ_stages.index(common_stage)
    CG4319_index = CG4319_stages.index(common_stage)
    CG12284_index = CG12284_stages.index(common_stage)
    values = [LacZ_values[lacZ_index], CG4319_values[CG4319_index], CG12284_values[CG12284_index]]
    samples = ["LacZ", "CG4319", "CG12284"]
    plt.bar(samples, values)
    plt.title("各系統の顕微鏡画像の輝度")
    plt.xlabel("系統")
    plt.ylabel("輝度(a.u.)")
    plt.savefig(output_path + "graph_stage" + str(common_stage) + ".jpg", dpi=300)
    plt.show()
