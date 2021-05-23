import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
data = pd.read_csv("Glucose_Desc.csv")
sns.set_context("paper", font_scale=1.8);
sns.set_style("ticks");
plt.figure(figsize=(6, 5.5));

reversed_data = data.iloc[::-1]

sns.scatterplot(x = "O/C", y = "H/C", palette="bright", marker="o", size="Generation", sizes=(90, 90), alpha=.90, hue = "Generation", legend="full",  data=reversed_data);
plt.legend(fontsize = 12, \
               #bbox_to_anchor= (1.03, 0.9), \
               title="Generation", \
               title_fontsize = 13, \
               facecolor = 'white');
# Set title
plt.title('')
# Set x-axis label
plt.xlabel('O/C');
# Set y-axis label
plt.ylabel('H/C');
#plt.xscale('log');
#plt.yscale('log');
plt.ylim(0,3.5);
plt.xlim(0,3.5);
plt.show()
grid = sns.FacetGrid(data, col = "Generation", col_wrap=3, palette="rocket", height=3, aspect=1);
grid.map(sns.scatterplot, "O/C", "H/C");
grid.set_axis_labels("O/C", "H/C");
grid.add_legend()


