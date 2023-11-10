import matplotlib.pyplot as plt
import japanize_matplotlib
import seaborn as sns


values = [25.775, 13.359, 19.545]
stages = ["LacZ", "CG4319", "CG12284"]
plt.bar(stages, values)
plt.title("各系統の顕微鏡画像の平均輝度")
plt.xlabel("系統")
plt.ylabel("平均輝度(a.u.)")
plt.savefig("/Users/kotaro/Desktop/graph.jpg")
plt.show()