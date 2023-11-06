
import pandas as pd
from scipy.stats import tukey_hsd
# import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# excel_path = "/Users/kotaro/Desktop/2023-10-31_実習B_7-12班.xls"

# df1 = pd.read_excel(excel_path,sheet_name="Sample Setup",skiprows=46)
# df2 = pd.read_excel(excel_path,sheet_name="Amplification Data",skiprows=46)
# df3 = pd.read_excel(excel_path,sheet_name="Results",skiprows=46)
# df4 = pd.read_excel(excel_path,sheet_name="Melt Curve Raw Data",skiprows=46)

# print(df1)
# print(df1.columns)

# print(df2)
# print(df2.columns)

# print(df3)
# print(df3.columns)

# print(df4)
# print(df4.columns)

test_path = "/Users/kotaro/Desktop/遺伝R/Test.csv"

df = pd.read_csv(test_path)

sns.violinplot(data=df, x="Condition", y="Value")#, hue="alive", fill=False)
plt.savefig("/Users/kotaro/Desktop/遺伝R/violin_test.jpg")
plt.show()
# plt.close()
sns.catplot(data=df, kind="swarm", x="Condition", y="Value")#, hue="smoker")
plt.savefig("/Users/kotaro/Desktop/遺伝R/swarm_test.jpg")
plt.show()
# plt.close()

print(df)

# group = df.groupby("Condition")
#
# temps = []
#
# for i,v in group:
#     # print(i)
#     temp = group.get_group(i)
#     temps.append(list(temp["Value"]))
#
# res = tukey_hsd(temps[0],temps[1],temps[2],temps[3])
# print(res.confidence_interval(confidence_level=0.99))
# print(res)


    