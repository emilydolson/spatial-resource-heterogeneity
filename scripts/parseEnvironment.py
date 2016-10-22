import sys, glob
from math import sqrt, log
from copy import deepcopy
import scipy.stats as stats

# Script to extract data about the spatial configuration of the environment from all of the environment files.
# Generates ents.csv, which is used by R scripts for stats.

# Most of this code is now part of the avidaspatial Python library.

def main():
    envs = parse_environment_files("environmentFiles/env*.cfg")
    outfile = open("../final_analysis/ents.csv", "w")
    outfile.write("seed, ent, functional_ent, overlap, variance, kurtosis, skew, resources\n")
    for key in envs.keys():
        print key
        data = envs[key]
        resList = [d[0] for d in data]
        rads = [d[1] for d in data]
        anchors = [d[2] for d in data]
        ent = calcEntropy(resList, anchors, rads, 59)
        functional_ent = calc_functional_entropy(resList, anchors, rads, 59)
        overlap = avg_overlap(resList, anchors, rads, 59)
        variance = overlap_variance(resList, anchors, rads, 59)
        kurtosis = overlap_kurtosis(resList, anchors, rads, 59)
        skew = overlap_skew(resList, anchors, rads, 59)
        niche_size = avg_niche_size(resList, anchors, rads, 59)
        resources = total_resources(resList, anchors, rads, 59)
        outfile.write(str(key) + ", " + str(ent)+"," + str(functional_ent) + "," + str(overlap) + "," + str(variance) + ", " + str(kurtosis) + ", " + str(skew) + "," + str(resources) +"\n")
    
    outfile.close()
         

def parse_environment_files(regex):
    names = glob.glob(regex)

    envs = {}
    for name in names:
        envs[name.strip("/environmentFls.cfg")] = parse_environment_file(name)

    return envs

def parse_environment_file(filename):
    infile = open(filename)
    lines = infile.readlines()
    infile.close()

    circles = []
    for line in lines:
        if line.startswith("GRADIENT_RESOURCE"):
            circles.append(parse_gradient(line))

    return circles
        

def parse_gradient(line):
    #remove "GRADIENT_RESOURCE"
    line = line[18:]

    sline = line.split(":")
    name = sline[0]
    radius = None
    x = None
    y = None

    for item in sline:
        if item.startswith("height"):
            radius = int(item.split("=")[1])
        elif item.startswith("peakx"):
            x = int(item.split("=")[1])
        elif item.startswith("peaky"):
            y = int(item.split("=")[1])

    return (name, radius, (x,y))


def makeNicheGrid(resList, anchors, rad, world_size):
    """
    Helper function for spatial heterogeneity calculations
    """
    resList = deepcopy(resList)
    #Reduce resource names to tasks being rewarded
    for i in range(len(resList)):
        resList[i] = resList[i][3:-1].lower()
        while resList[i][-1].isdigit():
            resList[i] = resList[i][:-1]

    if type(rad) == int:
        rad = [rad for i in range(len(anchors))]

    world = []
    for i in range(world_size):
        world.append([])
        for j in range(world_size):
            world[i].append(set())

    #Fill in data on niches present in each cell of the world
    for i in range(world_size):
        for j in range(world_size):
            for k in range(len(anchors)):
                if (dist((i,j), anchors[k])) <= rad[k]-1:
                    world[i][j].add(resList[k])
    return world

def overlap_skew(resList, anchors, rad, world_size, excludeDesert=True):
    niches = niche_analysis(resList, anchors, rad, world_size, excludeDesert)
    vals = []
    for key in niches.keys():
        vals.append(len(key)*(niches[key]/float(world_size*world_size)))

    return stats.skew(vals)

def overlap_kurtosis(resList, anchors, rad, world_size, excludeDesert=True):
    niches = niche_analysis(resList, anchors, rad, world_size, excludeDesert)
    vals = []
    for key in niches.keys():
        vals.append(len(key)*(niches[key]/float(world_size*world_size)))

    return stats.kurtosis(vals)

def overlap_variance(resList, anchors, rad, world_size, excludeDesert=True):
    niches = niche_analysis(resList, anchors, rad, world_size, excludeDesert)
    vals = []
    for key in niches.keys():
        vals.append(len(key)*(niches[key]/float(world_size*world_size)))

    return stats.tvar(vals)

def avg_overlap(resList, anchors, rad, world_size, excludeDesert=True):
    niches = niche_analysis(resList, anchors, rad, world_size, excludeDesert)
    total = 0.0
    for key in niches.keys():
        total += len(key)*(niches[key]/float(world_size*world_size))

    return total/len(niches.keys())

def avg_niche_size(resList, anchors, rad, world_size, excludeDesert=True):
    niches = niche_analysis(resList, anchors, rad, world_size, excludeDesert)
    total = 0.0
    for key in niches.keys():
        total += niches[key]

    return total/len(niches.keys())

def calc_functional_entropy(resList, anchors, rad, world_size, cutoff=36, excludeDesert=False):
    niches = niche_analysis(resList, anchors, rad, world_size, excludeDesert)
    for key in niches.keys():
        if niches[key] < cutoff:
            del niches[key]
    return entropy(niches)

def total_resources(resList, anchors, rad, world_size, excludeDesert=True):
    niches = niche_analysis(resList, anchors, rad, world_size, excludeDesert)
    total = 0.0
    for key in niches:
        total += len(key) * niches[key]

    return total


def calcEntropy(resList, anchors, rad, world_size, excludeDesert=False):
    return entropy(niche_analysis(resList, anchors, rad, world_size, excludeDesert))

def niche_analysis(resList, anchors, rad, world_size, excludeDesert=False):
    """
    Calculate the Shannon entropy of a given environment, treating each niche (where niches are defined by regions in which
    different sets of resources are rewarded) as a category. The environment is specified with the following inputs:
    
    resList - a list of strings of resource names. These names are expected to start with the letters "res" and end with numbers
    differentiating different copies of the same task. The middle characters denote the task being rewarded. For example, the first
    AND resource would be called resAND0.

    anchors - a list of tuples of ints denoting x,y anchor points in the world for each resource.

    radius - an int indicating the radius of the resource patches in the world.

    excludeDesert - an optional argument which defaults to False. If True is specific, niches in which no tasks are rewarded
    will not be considered in the calculation.
    """

    #Initialize list of list of sets to record which niches are where
    world = makeNicheGrid(resList, anchors, rad, world_size)
    
    #loop through world, counting frequency of each niche
    niches = {}
    for i in range(world_size):
        for j in range(world_size):
            if niches.has_key(frozenset(world[i][j])): #use frozensets because they are hashable
                niches[frozenset(world[i][j])] += 1
            else:
                niches[frozenset(world[i][j])] = 1

    if excludeDesert and niches.has_key(frozenset([])):
        del niches[frozenset([])]

    return niches

def dist(p1, p2):
    """
    Returns the distance between the two given tuples.
    """
    return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def entropy(dictionary):
    """
    Helper function for entropy calculations.
    Takes a frequency dictionary and calculates entropy of the keys.
    """
    total = 0.0
    entropy = 0
    for key in dictionary.keys():
        total += dictionary[key]

    for key in dictionary.keys():
        entropy += dictionary[key]/total * log(1.0/(dictionary[key]/total), 2)
    return entropy

if __name__ == "__main__":
    main()
