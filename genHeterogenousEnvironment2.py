#!/usr/bin/python
import argparse
import random
import sys, os
from math import sqrt, log, floor, ceil
import numpy as np
import scipy.stats as stats

#hpcc doesn't have matplotlib and ipythonblocks
from ipythonblocks import BlockGrid
import matplotlib
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Generate environment files for environmental heterogeneity experiment.", add_help=False)

parser.add_argument("--distance", default=None, type=int, help="Distance between anchor points of central patches")
parser.add_argument("--infinite", action="store_true", help="Use non-depletable resources")
parser.add_argument("--inflow", default=100, type=int, help="Sets amount of inflow.")
parser.add_argument("--outflow", default=.01, type=float, help="Sets amount of outflow.")
parser.add_argument("--patchesPerSide", default=4, type=int, help="number of patches in one row or column of world (assume square layout)")
parser.add_argument("--worldSize", default=59, type=int, help="Specifies size of one side (assumes square world)")
parser.add_argument("--backgroundOff", action="store_true", help="Don't put any resources in the desert")
parser.add_argument("--unlimitedEQU", action="store_true", help="Don't limit EQU (included for backwards compatability)")
parser.add_argument("--twoCircles", action="store_true", help="Triggers twocircles mode, which generates a set of environment files designed for the two circle overlaps experiments.")
parser.add_argument("--randAnchors", action="store_true", help="Generates random anchor points")

class args(object):
    pass

args = args()

parser.parse_args([],namespace=args)

frac = .0025 #fraction of available resource to use for reaction
resMax = 25 #maximum units of resource an organism can use per reaction
rxnType = "pow" #how is task reward granted?
maxCount = 1 #maximum number of times a task can be completed

gradResources = ["resNOT", "resAND", "resOR", "resNOR", "resNAND", "resORN", "resANDN", "resXOR"] #types of resources available in gradients

taskValDict = {"not":1.0,
               "nand":1.0,
               "and":2.0,
               "orn":2.0,
               "or":3.0,
               "andn":3.0,
               "nor":4.0,
               "xor":4.0,
               "equ":5.0} #rewards for each task



def main():

    mainParser = argparse.ArgumentParser(parents=[parser])

    mainParser.add_argument("outfile")
    mainParser.add_argument("randomSeed", default=0, type=int, help="The random seed to use.")
    mainParser.add_argument("patchRadius", default=5, type=int, help="Radius of resource patches")

    mainParser.parse_args(namespace=args)
    random.seed(args.randomSeed)

    background = True #is there always a background quantity of all resources
    if args.backgroundOff:
        background = False
  
    equLimited = True #is EQU limited
    if args.unlimitedEQU:
        equLimited = False

    if args.twoCircles:
        print "Generating", len(gradResources)*len(gradResources), "environment files."
        genTwoCircles()
        print "Success!"
        exit(0)
    
    outfile = open(args.outfile, "w")
  
    dist = (args.worldSize+1)/(args.patchesPerSide+1)
    
    anchors = []
    if args.randAnchors:
        anchors = calcRandomAnchors()
    elif args.distance == None:
        for i in range(dist-1, args.worldSize, dist):
            for j in range(dist-1, args.worldSize, dist):
                anchors.append((i,j))
    else:
        anchors = calcTightAnchors(args.distance, args.patchesPerSide)

    print anchors
    randResources = []
    nEach = args.patchesPerSide*args.patchesPerSide / len(gradResources)
    extras = args.patchesPerSide*args.patchesPerSide % len(gradResources)
    for i in range(nEach):
        for res in gradResources:
            randResources.append(res + str(i))

    additional = random.sample(gradResources, extras)
    for res in additional:
        randResources.append(res + str(nEach))

    random.shuffle(randResources)

    for i in range(len(anchors)):
        outfile.write(genGradient(randResources[i], args.inflow, args.patchRadius, anchors[i]))
    
    outfile.write("\n")
    outfile.write(genRes("resEQU1", args.inflow, args.outflow))
    outfile.write("\n")

    if background: #background resources
        for res in gradResources:
            outfile.write(genRes(res+"b", 1, .1))

        outfile.write("\n")

    for res in randResources:
        outfile.write(genReaction(res, not args.infinite))

    outfile.write(genReaction("resEQU1", equLimited))
    outfile.write("\n")

    if background: #background reactions
        for res in gradResources:
            outfile.write(genReaction(res+"b", 1))

        outfile.write("\n")
    ent = calcEntropy(randResources, anchors, args.patchRadius)

    outfile.write("# Entropy: " + str(ent))

    outfile.close()
    #graphEntropyVsDistance(25, randResources)
    #graphEntropyVsRadius(randResources, anchors)

def genTwoCircles():
    """
    Generates all 64 potential 2 circles comparison environments.
    """

    if args.distance == None:
        print "Error: genTwoCircles requires a distance"
        return
    
    centerPoint = (int(args.worldSize/2), int(args.worldSize/2))
    anchors = [(centerPoint[0]-int(floor(args.distance/2.0)), centerPoint[1]), (centerPoint[0]+int(ceil(args.distance/2.0)), centerPoint[1])]
    for res1 in gradResources:
        for res2 in gradResources:
            outfile = open(res1.strip("res").lower()+"X"+res2.strip("res").lower()+str(args.distance)+".cfg", "w")
            outfile.write(genGradient(res1+"1", args.inflow, args.patchRadius, anchors[0]))
            outfile.write(genGradient(res2+"2", args.inflow, args.patchRadius, anchors[1]))
            outfile.write("\n")
            
            outfile.write(genReaction(res1+"1", not args.infinite))
            outfile.write(genReaction(res2+"2", not args.infinite))

            outfile.close()            

def calcRandomAnchors():
    anchors = []
    rng = (args.patchRadius, args.worldSize - args.patchRadius)
    for i in range(args.patchesPerSide*args.patchesPerSide):
        anchors.append((random.randrange(rng[0], rng[1]), random.randrange(rng[0], rng[1])))

    return anchors

def calcTightAnchors(d, patches):
    centerPoint = (int(args.worldSize/2), int(args.worldSize/2))
    anchors = []
    if patches == 0:
        return anchors
    elif patches == 1:
        anchors.append(centerPoint)
        return anchors
   
    elif patches == 2:
        xs = [centerPoint[0]+int(floor(d/2)), centerPoint[0]-int(ceil(d/2))]
        ys = [centerPoint[1]+int(floor(d/2)), centerPoint[1]-int(ceil(d/2))]
        for i in xs:
            for j in ys:
                anchors.append((i,j))
        return anchors
    
    elif patches == 3:
        xs = [centerPoint[0] + int(d/2), centerPoint[0], centerPoint[0] - int(d/2)]
        ys = [centerPoint[1] + int(d/2), centerPoint[1], centerPoint[1] - int(d/2)]
        for i in xs:
            for j in ys:
                if not(i==centerPoint[0] and j==centerpoint[0]):
                    anchors.append((i,j))
        return anchors
    
    elif patches%2 == 0:
        dsout = (patches-2)/2 + 1
        xs = [centerPoint[0] + int(floor(d/2.0))+d*i for i in range(dsout)] + [centerPoint[0] - (int(ceil(d/2.0))+d*i) for i in range(dsout)]
        ys = [centerPoint[1] + int(floor(d/2.0))+d*i for i in range(dsout)] + [centerPoint[1] - (int(ceil(d/2.0))+d*i) for i in range(dsout)]
        for i in xs:
            anchors.append((i, max(ys)))
            anchors.append((i, min(ys)))
        for i in ys:
            anchors.append((max(xs), i))
            anchors.append((min(xs), i))

        if d != 0:
            anchors = list(set(anchors))
        anchors.sort()
        return (anchors + calcTightAnchors(d, patches-2))[:patches*patches] #to cut off the extras in the case where d=0

    else:
        #Note - an odd number of args.patchesPerSide requires that there be a patch at the centerpoint
        dsout = (patches-1)/2
        xs = [centerPoint[0] + d*i for i in range(dsout)] + [centerPoint[0] - d*i for i in range(dsout)]
        ys = [centerPoint[1] + d*i for i in range(dsout)] + [centerPoint[1] - d*i for i in range(dsout)]
        for i in xs:
            anchors.append((i, max(ys)))
            anchors.append((i, min(ys)))
        for i in ys:
            anchors.append((max(xs), i))
            anchors.append((min(xs), i))

        return anchors + calcTightAnchors(d, patches-2)


def genGradient(resource, inflow, radius, loc):
    return "GRADIENT_RESOURCE " + str(resource) + ":height=" + str(radius) + ":plateau=" + str(inflow) +":spread=" + str(radius-1) + ":updatestep=1000000:peakx=" + str(loc[0]) + ":peaky=" + str(loc[1]) + ":plateau_inflow=" + str(inflow) + ":initial=" + str(inflow) + "\n"

def genRes(resource, inflow, outflow):
    return "RESOURCE " + resource + ":inflow=" + str(inflow) + \
                ":outflow=" + str(outflow) + "\n"

def genReaction(resource, depletable=0):
    task = resource[3:-1].lower()
    name = resource[3:]
    return "REACTION " + name + " " + task + " process:resource=" + \
            resource + ":value=" + str(taskValDict[task]) + ":type=" \
            + rxnType + ":frac=" + str(frac) + ":max=" + str(resMax) + \
            ":depletable=" + str(int(depletable)) + " requisite:max_count=" \
            + str(maxCount) + "\n"

def graphEntropyVsDistanceVsRadius(resList):
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    plt.xlabel("Distance")
    plt.ylabel("Radius")
    plt.zlabel("Entropy")
    results = []
    errs = []
    for i in range(0, 60, 5):
        for j in range(0, 30, 5):
            ents = []
            anchors = calcTightAnchors(j, args.patchesPerSide)
            for k in range(30):
                random.shuffle(resList)
                ents.append(calcEntropy(resList, anchors, i))
            entAvg = sum(ents)/len(ents)
            results.append(entAvg)
            errs.append(np.sem(ents))

    plt.errorbar(results, yerr=errs)
    plt.savefig("entropyvsdistance.png")

def graphEntropyVsDistance(radius, resList):

    plt.xlabel("Distance")
    plt.ylabel("Entropy")
    results = []
    errs = []
    for i in range(60):
        ents = []
        anchors = calcTightAnchors(i, args.patchesPerSide)
        for j in range(2):
            random.shuffle(resList)
            ents.append(calcEntropy(resList, anchors, radius, False))
        entAvg = float(sum(ents))/len(ents)
        results.append(entAvg)
        errs.append(stats.sem(ents))

    plt.errorbar(range(len(results)), results, yerr=errs)
    plt.savefig("entropyvsdistance_radius" + str(radius)+ "_worldsize199.png")

def graphEntropyVsRadius(resList, anchors):

    plt.xlabel("Patch Radius")
    plt.ylabel("Entropy")
    results = []
    errs = []
    for i in range(60):
        ents = []
        for j in range(30):
            random.shuffle(resList)
            ents.append(calcEntropy(resList, anchors, i, True))
        entAvg = float(sum(ents))/len(ents)
        results.append(entAvg)
        errs.append(stats.sem(ents))

    plt.errorbar(range(len(results)), results, yerr=errs)
    plt.savefig("entropyvsradius_errors_nodesert.png")

def calcEntropy(resList, anchors, rad, excludeDesert=False):
    world = []
    for i in range(args.worldSize):
        world.append([])
        for j in range(args.worldSize):
            world[i].append(set())

    for i in range(args.worldSize):
        for j in range(args.worldSize):
            for k in range(len(anchors)):
                if (dist((i,j), anchors[k])) <= rad-1:
                    world[i][j].add(resList[k][3:-1])
    
    entropy = 0
    niches = {}
    for i in range(args.worldSize):
        for j in range(args.worldSize):
            if niches.has_key(frozenset(world[i][j])):
                niches[frozenset(world[i][j])] += 1
            else:
                niches[frozenset(world[i][j])] = 1

    if excludeDesert and niches.has_key(frozenset([])):
        del niches[frozenset([])]

    total = 0.0
    for key in niches.keys():
        total += niches[key]

    for key in niches.keys():
        entropy += niches[key]/total * log(1.0/(niches[key]/total), 2)
    
    return entropy

def nTasks(decNum):
    bitstring = bin(int(decNum))[2:]
    return bitstring.count("1")

def makeNTasksImage(filename):
    data = loadGridData(filename, "nTasks")
    grid = []

    for row in range(len(data)):
        grid.append([])
        for col in range(len(data[row])):
            arr = np.zeros((1,1,3))
            if float(data[row][col]) != 0:
                arr[0,0,0] = (float(data[row][col])/9.0)
            else:
                arr[0,0,0] = 0

            arr[0,0,1] = 1
            arr[0,0,2] = 1
            rgb = matplotlib.colors.hsv_to_rgb(arr)
            grid[row].append(rgb)
    plt.imshow(grid)
            
    return grid

def loadGridData(filename, mode):
    infile = open(filename)
    data = []
    for line in infile:
        if mode == "nTasks":
            data.append([nTasks(s) for s in line.split()])
        elif mode == "species":
            data.append([((log(float(s)+1, 2)/11)) for s in line.split()])
    infile.close()
    return data

def loadAggregateGridData(filelist, mode):
    data = []
    for i in range(args.worldSize):
        data.append([])
        for j in range(args.worldSize):
            if mode.startswith("mostcommon"):
                data[i].append([])
            else:
                data[i].append(0.0)

    for f in filelist:
        infile = open(f)
        lines = infile.readlines()
        for i in range(args.worldSize):
            if mode.find("nTasks") != -1:
                lines[i] = [nTasks(num)/9.0 for num in lines[i].split()]
            elif mode.find("species") != -1:
                lines[i] = [log(float(num)+1, 2)/11.0 for num in lines[i].split()]
            else:
                print "Warning: using raw numbers"
                lines[i] = lines[i].split()
            for j in range(args.worldSize):
                if mode.startswith("mostcommon"):
                    data[i][j].append(float(lines[i][j]))
                else:
                    data[i][j] += float(lines[i][j])
        infile.close()

    if mode.startswith("mostcommon"):
        for i in range(args.worldSize):
            for j in range(args.worldSize):
                data[i][j] = float(max(set(data[i][j]), key=data[i][j].count))
    else:
        for i in range(args.worldSize):
            for j in range(args.worldSize):
                data[i][j] /= len(filelist)

    if mode.find("intensify") != -1:
        for i in range(args.worldSize):
            for j in range(args.worldSize):
                data[i][j] *= 3.0

    return data

def colorGrid(data):
    grid = BlockGrid(args.worldSize, args.worldSize)
    for row in range(len(data)):
        for col in range(len(data[row])):
            arr = np.zeros((1,1,3))
            if float(data[row][col]) != 0:
                arr[0,0,0] = float(data[row][col])
            else:
                arr[0,0,0] = 0

            arr[0,0,1] = 1
            arr[0,0,2] = 1

            rgb = matplotlib.colors.hsv_to_rgb(arr)
            data[row][col] = [rgb[0,0,0], rgb[0,0,1], rgb[0,0,2]]
            grid[row,col].red = rgb[0,0,0]*255
            grid[row,col].green = rgb[0,0,1]*255
            grid[row,col].blue = rgb[0,0,2]*255

    return grid, data

def makeMostCommonGrid(filelist, mode="species"):
    data = loadAggregateGridData(filelist, "mostcommon-"+mode)
    return colorGrid(data)

def makeNTasksGrid(filename):
    data = loadGridData(filename, "nTasks")
    return colorGrid(data)

def makeAggregateNTasksGrid(filelist, intensify=False):
    if intensify:
        data = loadAggregateGridData(filelist, "nTasks-intensify")
    else:
        data = loadAggregateGridData(filelist, "nTasks")
    return colorGrid(data)

def makeSpeciesGrid(filename):
    data = loadGridData(filename, "species")
    return colorGrid(data)

def makeAggregateSpeciesGrid(filelist):
    data = loadAggregateGridData(filelist, "species")
    return colorGrid(data)
    
def calcGrid(resList, anchors, rad):
    world = []
    grid = BlockGrid(args.worldSize, args.worldSize)
    for i in range(args.worldSize):
        world.append([])
        for j in range(args.worldSize):
            world[i].append(set())

    for i in range(args.worldSize):
        for j in range(args.worldSize):
            for k in range(len(anchors)):
                if (dist((i,j), anchors[k])) <= rad-1:
                    world[i][j].add(resList[k][3:-1])
    
    entropy = 0
    niches = {}
    for i in range(args.worldSize):
        for j in range(args.worldSize):
            greenVal = 0.0
            blueVal = 0.0
            redVal = 0.0
            colorDict = {"AND":(255, 0, 0), "OR":(0, 255, 0), "ORN":(0,0,255), "ANDN":(255, 255, 0), "XOR":(255,0,255), "NOR":(0, 255, 255), "NOT":(255, 155, 0), "NAND":(255, 255, 255)}
            for res in world[i][j]:
                redVal += colorDict[res][0]
                greenVal += colorDict[res][1]
                blueVal += colorDict[res][2]
            if len(world[i][j]) > 0:
                redVal/=len(world[i][j])
                greenVal /= len(world[i][j])
                blueVal /= len(world[i][j])
            grid[i,j].red = redVal
            grid[i,j].green = greenVal
            grid[i,j].blue = blueVal
    grid.lines_on = False
    return grid

def calcNResGrid(resList, anchors, rad):
    world = []
    grid = BlockGrid(args.worldSize, args.worldSize)
    for i in range(args.worldSize):
        world.append([])
        for j in range(args.worldSize):
            world[i].append(0)

    for i in range(args.worldSize):
        for j in range(args.worldSize):
            for k in range(len(anchors)):
                if (dist((i,j), anchors[k])) <= rad-1:
                    world[i][j] += 1
    
    entropy = 0
    niches = {}
    for i in range(args.worldSize):
        for j in range(args.worldSize):
            arr = np.zeros((1,1,3))
            if float(world[i][j]) != 0:
                arr[0,0,0] = (float(world[i][j]/10.0)*.6)
            else:
                arr[0,0,0] = 0
                
            arr[0,0,1] = 1
            arr[0,0,2] = 1
            rgb = matplotlib.colors.hsv_to_rgb(arr)
            grid[i,j].red = rgb[0,0,0]*255
            grid[i,j].green = rgb[0,0,1]*255
            grid[i,j].blue = rgb[0,0,2]*255

    grid.lines_on = False
    return grid

def dist(p1, p2):
    return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

if __name__ == "__main__":
    main()
