import sys
import avidaspatial as avsp
import glob

d14 = glob.glob("/media/emily/hdd/resource-heterogeneity/experiment/"+sys.argv[1]+"/*/grid_task.100000.dat")
print(d14)
found = [0]*30
grid_14 = avsp.load_grid_data(d14)
for y in range(59):
    for x in range(59):
        for rep in range(30):
            if grid_14[y][x][rep] == "0b11":
                found[rep] = 1

print(sum(found))
