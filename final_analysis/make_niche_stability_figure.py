import avidaspatial as avsp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob
import seaborn as sns

d0 = glob.glob("/media/emily/hdd/resource-heterogeneity/experiment/2circles_notXorn0/*/grid_task.100000.dat")
d14 = glob.glob("/media/emily/hdd/resource-heterogeneity/experiment/2circles_notXorn14/*/grid_task.100000.dat")
d30 = glob.glob("/media/emily/hdd/resource-heterogeneity/experiment/2circles_notXorn30/*/grid_task.100000.dat")

grid_0 = avsp.agg_grid(avsp.load_grid_data(d0))
grid_14 = avsp.agg_grid(avsp.load_grid_data(d14))
grid_30 = avsp.agg_grid(avsp.load_grid_data(d30))

#print(grid_14)
ranks = {'0b11': 3, '0b10': 2, '0b1': 1, '0b0': 0, '0b00': 0, '0b01': 1}
#ranks = avsp.generate_ranks(grid_14, 6)
#print(ranks)

grid_0 = avsp.assign_ranks_to_grid(grid_0, ranks)
grid_14 = avsp.assign_ranks_to_grid(grid_14, ranks)
grid_30 = avsp.assign_ranks_to_grid(grid_30, ranks)

#print(grid_0)

pal = sns.cubehelix_palette(4, start=.5, rot=-.75)
#pal = sns.color_palette("cubehelix", 3)
pal.pop()
pal.insert(0, (0,0,0))

fig = plt.figure(figsize=(7.5, 2.5))

ax1 = fig.add_subplot(131)
avsp.make_imshow_plot(avsp.color_grid(grid_0, pal, 3, False), "teststability.png")
ax1.text(2, 7, "a", fontsize=20, color="white")
patches = []
labels = ["None", "NOT", "ORN", "Both"]
for i in range(len(pal)):
    patches.append(mpatches.Patch(facecolor=pal[i], edgecolor="white", linewidth=1, label = labels[i]))

l = plt.legend(handles=patches, ncol=2, loc="lower center")
for text in l.get_texts():
    text.set_color("white")
    text.set_size(14)

center = (29, 29)
r = 14

circle = plt.Circle(center, r, edgecolor="blue", linewidth=3, fc="none", zorder=3)
ax1.add_patch(circle)

#print ax1.patches

ax2 = fig.add_subplot(132)
avsp.make_imshow_plot(avsp.color_grid(grid_14, pal, 3, False), "teststability.png")
ax2.text(2, 7, "b", fontsize=20, color="white")

center1 = (29, 22)
circle1 = plt.Circle(center1, r, edgecolor="gold", linewidth=3, fc="none", zorder=3)
center2 = (29, 36)
circle2 = plt.Circle(center2, r, edgecolor="green", linewidth=3, fc="none", zorder=3)
ax2.add_patch(circle1)
ax2.add_patch(circle2)


ax3 = fig.add_subplot(133)
avsp.make_imshow_plot(avsp.color_grid(grid_30, pal, 3, False), "teststability.png")
ax3.text(2, 7, "c", fontsize=20, color="white")
center1 = (29, 14)
circle1 = plt.Circle(center1, r, edgecolor="gold", linewidth=3, fc="none", zorder=3)
center2 = (29, 44)
circle2 = plt.Circle(center2, r, edgecolor="green", linewidth=3, fc="none", zorder=3)
ax3.add_patch(circle1)
ax3.add_patch(circle2)

#plt.show()
plt.savefig("nicheConstruction.eps")

