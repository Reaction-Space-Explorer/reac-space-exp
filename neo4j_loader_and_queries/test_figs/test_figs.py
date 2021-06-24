import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(r"D:\BMSIS_YSP\reac-space-exp\neo4j_loader_and_queries\output\2021-01-26_17-18-42-588599\all_generations_abundance_scores.csv")



print(df.head())

df.plot.area(x='rank_node_deg_by_gen',
             y='count_relationships')

plt.show()