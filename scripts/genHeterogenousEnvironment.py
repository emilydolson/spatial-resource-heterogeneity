#!/usr/bin/python


import random
import sys, os
from math import sqrt, log, floor, ceil
import numpy as np
import scipy.stats as stats

#hpcc doesn't have matplotlib and ipythonblocks
#from ipythonblocks import BlockGrid
#import matplotlib
#from matplotlib import pyplot as plt

patchRadius = 5
patchesPerSide = 4 #number of patches in one row or column of world (assume square layout)
worldSize = 59 #Specifies size of one side (assumes square world)

background = True #is there always a background quantity of all resources
equLimited = True #is EQU limited

infinite = False #Are resources unlimited?
inflow = 100 #Resource inflow rate
outflow = .01 #Resource outflow rate
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
    randSeed = None
    distance = None
    global patchRadius

    if len(sys.argv) < 2:
        print "Usage: ./genHeterogenousEnvironment outputFileName [seed] [patchRadius] [distance]"
        return

    if len(sys.argv) == 3:
        randSeed = int(sys.argv[2])
    elif len(sys.argv) == 4:
        patchRadius = int(sys.argv[3])
        randSeed = int(sys.argv[2])
    elif len(sys.argv) == 5:
        patchRadius = int(sys.argv[3])
        randSeed = int(sys.argv[2])
        distance = int(sys.argv[4])

    random.seed(randSeed)

    outfile = open(sys.argv[1], "w")
  
    dist = (worldSize+1)/(patchesPerSide+1)
    
    anchors = []
    if distance == None:
        for i in range(dist-1, worldSize, dist):
            for j in range(dist-1, worldSize, dist):
                anchors.append((i,j))

    else:
        anchors = calcTightAnchors(distance, patchesPerSide)

    print anchors
    randResources = []
    nEach = patchesPerSide*patchesPerSide / len(gradResources)
    extras = patchesPerSide*patchesPerSide % len(gradResources)
    for i in range(nEach):
        for res in gradResources:
            randResources.append(res + str(i))

    additional = random.sample(gradResources, extras)
    for res in additional:
        randResources.append(res + str(nEach))

    random.shuffle(randResources)

    for i in range(len(anchors)):
        outfile.write(genGradient(randResources[i], inflow, patchRadius, anchors[i]))
    
    outfile.write("\n")
    outfile.write(genRes("resEQU1", inflow, outflow))
    outfile.write("\n")

    if background: #background resources
        for res in gradResources:
            outfile.write(genRes(res+"b", 1, .1))

        outfile.write("\n")

    for res in randResources:
        outfile.write(genReaction(res, not infinite))

    outfile.write(genReaction("resEQU1", equLimited))
    outfile.write("\n")

    if background: #background reactions
        for res in gradResources:
            outfile.write(genReaction(res+"b", 1))

        outfile.write("\n")
    ent = calcEntropy(randResources, anchors, patchRadius)

    outfile.write("# Entropy: " + str(ent))

    outfile.close()
    #graphEntropyVsDistance(25, randResources)
    #graphEntropyVsRadius(randResources, anchors)

def calcTightAnchors(d, patches):
    centerPoint = (int(worldSize/2), int(worldSize/2))
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
        #Note - an odd number of patchesPerSide requires that there be a patch at the centerpoint
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
            anchors = calcTightAnchors(j, patchesPerSide)
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
        anchors = calcTightAnchors(i, patchesPerSide)
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
    for i in range(worldSize):
        world.append([])
        for j in range(worldSize):
            world[i].append(set())

    for i in range(worldSize):
        for j in range(worldSize):
            for k in range(len(anchors)):
                if (dist((i,j), anchors[k])) <= rad-1:
                    world[i][j].add(resList[k][3:-1])
    
    entropy = 0
    niches = {}
    for i in range(worldSize):
        for j in range(worldSize):
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

def calcGrid(resList, anchors, rad):
    world = []
    grid = BlockGrid(worldSize, worldSize)
    for i in range(worldSize):
        world.append([])
        for j in range(worldSize):
            world[i].append(set())

    for i in range(worldSize):
        for j in range(worldSize):
            for k in range(len(anchors)):
                if (dist((i,j), anchors[k])) <= rad-1:
                    world[i][j].add(resList[k][3:-1])
    
    entropy = 0
    niches = {}
    for i in range(worldSize):
        for j in range(worldSize):
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

def dist(p1, p2):
    return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

if __name__ == "__main__":
    main()
