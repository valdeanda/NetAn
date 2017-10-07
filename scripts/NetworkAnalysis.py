#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     NetworkAnalysis.py
# Purpose:  Compute topological statistics and some figures from adyacence list  
# @uthor:   Marcos Emmanuel Gonzales Laffitte  - laffitte6345@live.com.m
# modifications: Valerie de Anda -vdeanda@ciencias.unam.mx
## Created:    September  2017
# ------------------------------                                                                                                                                                                                                                                                  #
import os
import argparse
from sys import argv
import networkx as nx
from numpy import arange
import math
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter("ignore")
plt.switch_backend('agg')
import community

#Options 
epilog = """Example:

$  python3 scripts/NetworkAnalysis.py -d data/a_phylum_consensus.txt """ 
          
parser = argparse.ArgumentParser(description=__doc__, epilog=epilog)
parser.add_argument('filename',
                    help="Input file in tabular format of adyacence list ")
parser.add_argument(
    '-u', '--undirect', help='Undirected graphs no weigth')

parser.add_argument('-d', '--direct', help='Directed graph') 
args= parser.parse_args() 

originalUndirNetwork = nx.Graph()
originalDirecNetwork = nx.DiGraph()
networkCounter = 0
networkType = ""

# Function: analize network ------------------------------
def CleanNetworkFilename(networkFileName):
    # variables
    nameArray = []
    name = ""
    if("/" in networkFileName):
        nameArray = networkFileName.split("/")
        name = nameArray[-1]
        return(name)
    if("\\" in networkFileName):
        nameArray = networkFileName.split("\\")
        name = nameArray[-1]
        return(name)
    return(networkFileName)

# Function: ParseNetwork ----------------------------------
def ParseFileToNetwork(networkFileName, typeNW, weight):
    # function message
    if(weight != "community"):
        print("\t- Parsing file to network ...")
    # variables
    FILE = None
    eachLine = None
    possEdge = []
    firstNode = ""
    secondNode = ""
    edgeWeight = 0
    # define a directed or undirected network
    if(typeNW == "-u"):
        fileNetwork = nx.Graph()
    if(typeNW == "-d"):
        fileNetwork = nx.DiGraph()
    # openfile
    FILE = open(networkFileName, "r")
    # create network looping to get the arrows
    for eachLine in FILE:
        possEdge = eachLine.split("\t")
        if(len(possEdge) >= 3):
            firstNode = possEdge[0].strip()
            secondNode = possEdge[1].strip()
            if(weight == "community"):
                edgeWeight = 1
            else:
                edgeWeight = float(possEdge[2].strip())
            # add edge to network
            fileNetwork.add_edge(firstNode, secondNode, weight = edgeWeight)
    FILE.close()
    # end of function
    return(fileNetwork)

# Function: analize order ---------------------------------
def AnalizeOrder(someNetwork):
    # function message
    print("\t- Obtaining order ...")
    # variables
    orderFinalResult = ""
    orderResult = "Order:"
    # get order
    orderResult = orderResult + str(someNetwork.order())
    # end of function
    orderFinalResult = "\n\n" + orderResult + "\n\n"
    return(orderFinalResult)

# Function: analize size -----------------------------------
def AnalizeSize(someNetwork):
    # function message
    print("\t- Obtaining size ...")
    # variables
    sizeFinalResult = ""
    sizeResult = "Size:"
    # get size
    sizeResult = sizeResult + str(someNetwork.size())
    # end of function
    sizeFinalResult = "\n\n" + sizeResult + "\n\n"
    return(sizeFinalResult)

# Function: analize diameter ----------------------------------
def AnalizeDiameter(someNetwork):
    # function message
    print("\t- Obtaining diameter (for undir network version) ...")
    # variables
    diameterFinalResult = ""
    diameterResult = "Diameter (Undirected):"
    # get diameter
    if(nx.is_connected(someNetwork.to_undirected())):
        diameterResult = diameterResult + str(nx.diameter(someNetwork.to_undirected()))
    else:
        diameterResult = diameterResult + "NOT CONNECTED"
    # end of function
    diameterFinalResult = "\n\n" + diameterResult + "\n\n"
    return(diameterFinalResult)

# Function: analize radius ------------------------------------
def AnalizeRadius(someNetwork):
    # function message
    print("\t- Obtaining radius (for undir network version) ...")
    # variables
    radiusFinalResult = ""
    radiusResult = "Radius (Undirected):"
    # get radius
    if(nx.is_connected(someNetwork.to_undirected())):
        radiusResult = radiusResult + str(nx.radius(someNetwork.to_undirected()))
    else:
        radiusResult = radiusResult + "NOT CONNECTED"
    # end of function
    radiusFinalResult = "\n\n" + radiusResult + "\n\n"
    return(radiusFinalResult)

# Function: analize density -------------------------------------
def AnalizeDensity(someNetwork):
    # function message
    print("\t- Obtaining density ...")
    # variables
    densityFinalResult = ""
    densityResult = "Density:"
    # get density
    densityResult = densityResult + str(nx.density(someNetwork))
    # end of function
    densityFinalResult = "\n\n" + densityResult + "\n\n"
    return(densityFinalResult)

# Function: analize mean degree ----------------------------------
def AnalizeMeanDegree(someNetwork, typeNW):
    # function message
    print("\t- Obtaining mean degree ...")
    # variables
    meanDegreeFinalResult = ""
    meanDegreeResult = "Mean Degree:"
    # get mean degree
    if(typeNW == "-u"):
        meanDegreeResult = meanDegreeResult + str((2*someNetwork.size())/(someNetwork.order()))
    if(typeNW == "-d"):
        meanDegreeResult = meanDegreeResult + str((someNetwork.size())/(someNetwork.order()))
    # end of function
    meanDegreeFinalResult = "\n\n" + meanDegreeResult + "\n\n"
    return(meanDegreeFinalResult)

# Function: analize max degree -----------------------------------
def AnalizeMaxDegree(someNetwork):
    # function message
    print("\t- Obtaining max degree ...")
    # variables
    maxDegreeFinalResult = ""
    maxDegreeResult = "Max Possible Degree:"
    # get max degree
    maxDegreeResult = maxDegreeResult + str((someNetwork.order()) - 1)
    # end of function
    maxDegreeFinalResult = "\n\n" + maxDegreeResult + "\n\n"
    return(maxDegreeFinalResult)

# Function: analize mean clustering coefficient -------------------
def AnalizeClustCoeff(someNetwork, typeNW):
    # function message
    print("\t- Obtaining mean clusttering coefficient ...")
    # variables
    clustCoeffFinalResult = ""
    clustCoeffResult = "Clustering Coefficient:"
    # get clustering coefficient
    if(nx.is_connected(someNetwork.to_undirected())):
        if(typeNW == "-u"):
            clustCoeffResult = clustCoeffResult + str(nx.average_clustering(someNetwork.to_undirected()))
        if(typeNW == "-d"):
            clustCoeffResult = clustCoeffResult + str((nx.average_clustering(someNetwork.to_undirected()))/2)
    else:
        clustCoeffResult = clustCoeffResult + "NOT CONNECTED"
    # end of function
    clustCoeffFinalResult = "\n\n" + clustCoeffResult + "\n\n"
    return(clustCoeffFinalResult)

# Function: analize max degree hubs ---------------------------------
def AnalizeMaxHubs(someNetwork, typeNW):
    # function message
    print("\t- Obtaining hubs ...")
    # variables
    maxHubsFinalResult = ""
    maxHubsResult = "Max Degree Hubs:\n"
    degMapping = dict()
    inDegMapping = dict()
    outDegMapping = dict()
    nodeVector = []
    inNodeVector = []
    outNodeVector = []
    degVector = []
    inDegVector = []
    outDegVector = []
    hubsVector = []
    inHubsVector = []
    outHubsVector = []
    maxDegree = 0
    inMaxDegree = 0
    outMaxDegree = 0
    counter = 0
    # get max degree hubs undirected graph
    if(typeNW == "-u"):
        degMapping = dict(someNetwork.degree())
        nodeVector = list(degMapping.keys())
        degVector = list(degMapping.values())
        degVector.sort()
        degVector.reverse()
        maxDegree = degVector[0]
        for counter in range(len(nodeVector)):
            if(degMapping[nodeVector[counter]] == maxDegree):
                hubsVector.append(nodeVector[counter])
        # get result string
        maxHubsResult = maxHubsResult + "- Max Degree: " + str(maxDegree) + "\n"
        maxHubsResult = maxHubsResult + "- Number of Hubs with Max Degree: " + str(len(hubsVector)) + "\n"
        maxHubsResult = maxHubsResult + "- Hubs with Max Degree: " + ",".join(hubsVector)
    # get max degree hubs directed graph
    if(typeNW == "-d"):
        inDegMapping = dict(someNetwork.in_degree())
        outDegMapping = dict(someNetwork.out_degree())
        inNodeVector = list(inDegMapping.keys())
        outNodeVector = list(outDegMapping.keys())
        inDegVector = list(inDegMapping.values())
        outDegVector = list(outDegMapping.values())
        inDegVector.sort()
        inDegVector.reverse()
        outDegVector.sort()
        outDegVector.reverse()
        inMaxDegree = inDegVector[0]
        outMaxDegree = outDegVector[0]
        for counter in range(len(inNodeVector)):
            if(inDegMapping[inNodeVector[counter]] == inMaxDegree):
                inHubsVector.append(inNodeVector[counter])
        for counter in range(len(outNodeVector)):
            if(outDegMapping[outNodeVector[counter]] == outMaxDegree):
                outHubsVector.append(outNodeVector[counter])
        # get result string
        maxHubsResult = maxHubsResult + "- Max In Degree: " + str(inMaxDegree) + "\n"
        maxHubsResult = maxHubsResult + "- Number of Hubs with Max In Degree: " + str(len(inHubsVector)) + "\n"
        maxHubsResult = maxHubsResult + "- Hubs with Max In Degree: " + ",".join(inHubsVector) + "\n"
        maxHubsResult = maxHubsResult + "- Max Out Degree: " + str(outMaxDegree) + "\n"
        maxHubsResult = maxHubsResult + "- Number of Hubs with Max Out Degree: " + str(len(outHubsVector)) + "\n"
        maxHubsResult = maxHubsResult + "- Hubs with Max Out Degree: " + ",".join(outHubsVector)
    # end of function
    maxHubsFinalResult = "\n\n" + maxHubsResult + "\n\n"
    return(maxHubsFinalResult)

# Function: analize minimum cut vertices -------------------------
def AnalizeCutVertices(someNetwork):
    # function message
    print("\t- Obtaining mincut vertex set (for undir network version) ...")
    # variables
    cutVerticesFinalResult = ""
    cutVerticesResult = "Minimum Cut Node Set (Undirected):\n"
    cutVerticesVector = []
    # get cut vertices
    if(nx.is_connected(someNetwork.to_undirected())):
        cutVerticesVector = list(nx.minimum_node_cut(someNetwork.to_undirected()))
        cutVerticesResult = cutVerticesResult + "- Number of Vertices in MinCut Set: " + str(len(cutVerticesVector)) + "\n"
        cutVerticesResult = cutVerticesResult + "- Vertices in MinCut Set: " + ",".join(cutVerticesVector)
    else:
        cutVerticesResult = cutVerticesResult + "NOT CONNECTED"
    # end of function
    cutVerticesFinalResult = "\n\n" + cutVerticesResult + "\n\n"
    return(cutVerticesFinalResult)

# Function: analize cut edges -------------------------------------
def AnalizeCutEdges(someNetwork):
    # function message
    print("\t- Obtaining mincut edge set (for undir network version) ...")
    # variables
    cutEdgesFinalResult = ""
    cutEdgesResult = "Minimum Cut Edge Set (Undirected):\n"
    cutEdgesVector = []
    cutEdgesStr = ""
    # get cut edges
    if(nx.is_connected(someNetwork.to_undirected())):
        cutEdgesVector = list(nx.minimum_edge_cut(someNetwork.to_undirected()))
        cutEdgesStr = str(cutEdgesVector).strip("[]")
        cutEdgesResult = cutEdgesResult + "- Number of Edges in MinCut Set: " + str(len(cutEdgesVector)) + "\n"
        cutEdgesResult = cutEdgesResult + "- Edges in MinCut Set: " + cutEdgesStr
    else:
        cutEdgesResult = cutEdgesResult + "NOT CONNECTED"
    # end of function
    cutEdgesFinalResult = "\n\n" + cutEdgesResult + "\n\n"
    return(cutEdgesFinalResult)

# Function: analize connected components ------------------------------
def AnalizeConnComps(someNetwork):
    # function message
    print("\t- Obtaining connected components ...")
    # variables
    connCompsFinalResult = ""
    connCompsResult = "Connected Components:\n"
    connCompsVector = []
    numbConnComps = 0
    counter = 0
    # get connected components
    connCompsVector = list(nx.connected_components(someNetwork.to_undirected()))
    numbConnComps = len(connCompsVector)
    # get result string
    connCompsResult = connCompsResult + "- Number of CC's: " + str(numbConnComps) + "\n"
    for counter in range(numbConnComps):
        connCompsResult = connCompsResult + "- Number of nodes in CC_" + str(counter + 1) + " : " + str(len(list(connCompsVector[counter]))) + "\n"
        connCompsResult = connCompsResult + "- Nodes in CC_" + str(counter + 1) + ": " + ",".join(list(connCompsVector[counter])) + "\n"
    # end of function
    connCompsFinalResult = "\n\n" + connCompsResult + "\n\n"
    return(connCompsFinalResult)

# Function: analize strongly connected components ----------------------
def AnalizeStronConnComps(someNetwork, typeNW):
    # function message
    print("\t- Obtaining strongly connected components (just for directed case) ...")
    # variables
    stronConnCompsFinalResult = ""
    stronConnCompsResult = "Strongly Connected Components:\n"
    stronConnCompsVector = []
    numbStronConnComps = 0
    counter = 0
    # get strongly connected components
    if(nx.is_directed(someNetwork)):
        stronConnCompsVector = list(nx.strongly_connected_components(someNetwork))
        numbStronConnComps = len(stronConnCompsVector)
        # get result string
        stronConnCompsResult = stronConnCompsResult + "- Number of SCC's: " + str(numbStronConnComps) + "\n"
        for counter in range(numbStronConnComps):
            stronConnCompsResult=stronConnCompsResult+"- Number of nodes in SCC_"+str(counter+1)+" : "+str(len(list(stronConnCompsVector[counter])))+"\n"
            stronConnCompsResult = stronConnCompsResult + "- Nodes in SCC_" + str(counter + 1) + ": " + ",".join(list(stronConnCompsVector[counter])) + "\n"
    else:
        stronConnCompsResult = stronConnCompsResult + "NOT DIRECTED"
    # end of function
    stronConnCompsFinalResult = "\n\n" + stronConnCompsResult + "\n\n"
    return(stronConnCompsFinalResult)

# Function: analize max cliques --------------------------------------
def AnalizeMaxCliques(someNetwork):
    # function message
    print("\t- Obtaining maximal cliques (for undir network version) ...")
    # variables
    maxCliquesFinalResult = ""
    maxCliquesResult = "Max Cliques:\n"
    maxCliquesVector = []
    numbMaxCliques = 0
    counter = 0
    # get max cliques
    maxCliquesVector = list(nx.find_cliques(someNetwork.to_undirected()))
    numbMaxCliques = len(maxCliquesVector)
    # get result string
    maxCliquesResult = maxCliquesResult + "- Number of Maximal Cliques: " + str(numbMaxCliques) + "\n"
    for counter in range(numbMaxCliques):
        maxCliquesResult = maxCliquesResult + "- Number of nodes in Maximal Clique " + str(counter + 1) + " : " + str(len(maxCliquesVector[counter])) + "\n"
        maxCliquesResult = maxCliquesResult + "- Nodes in Clique_" + str(counter + 1) + ": " + ",".join(maxCliquesVector[counter]) + "\n"
    # end of function
    maxCliquesFinalResult = "\n\n" + maxCliquesResult + "\n\n"
    return(maxCliquesFinalResult)

# Function: analize cycle basis --------------------------------------
def AnalizeCycleBasis(someNetwork):
    # function message
    print("\t- Obtaining cycle basis (for undir network version) ,,,")
    # variables
    basicCyclesFinalResult = ""
    basicCyclesResult = "Cycle Basis (Undirected):\n"
    basicCyclesVector = []
    numbBasicCycles = 0
    counter = 0
    # get cycle basis
    basicCyclesVector = list(nx.cycle_basis(someNetwork.to_undirected()))
    numbBasicCycles = len(basicCyclesVector)
    # get result string
    basicCyclesResult = basicCyclesResult + "- Number Cycles in CB: " + str(numbBasicCycles) + "\n"
    for counter in range(numbBasicCycles):
        basicCyclesResult = basicCyclesResult + "- Number of nodes in Cycle " + str(counter + 1) + ": " + str(len(basicCyclesVector[counter])) + "\n"
        basicCyclesResult = basicCyclesResult + "- Nodes in Cycle_" + str(counter + 1) + ": " + ",".join(basicCyclesVector[counter]) + "\n"
    # end of function
    basicCyclesFinalResult = "\n\n" + basicCyclesResult + "\n\n"
    return(basicCyclesFinalResult)

# Function: analize independent set -----------------------------------
def AnalizeIndependentSets(someNetwork, totIndSet):
    # function message
    print("\t- Obtaining maximal independent set ...")
    # variablesposs
    independentSetsFinalResult = ""
    independentSetsResult = "Maximal Independet Set (might be repeated):\n"
    maxIndependentSet = []
    counter = 0
    # get independent sets
    for counter in range(totIndSet):
        independentSetsResult = independentSetsResult + "- Independent Set: " + str(counter + 1) + "\n"
        maxIndependentSet = nx.maximal_independent_set(someNetwork.to_undirected())
        independentSetsResult = independentSetsResult + "- Number of nodes in this maximal independent set: " + str(len(maxIndependentSet)) + "\n"
        independentSetsResult = independentSetsResult + "- Nodes in this maximal independent set: " + ",".join(maxIndependentSet) + "\n"
    # end of function
    independentSetsFinalResult = "\n\n" + independentSetsResult + "\n\n"
    return(independentSetsFinalResult)

# Function: analize communities -----------------------------------------
def AnalizeCommunitiesAndDrawings(fileName, someNetwork):
    # function message
    print("\t- Obtaining communities (for undir network version)...")
    # variables
    communitiesFinalResult = ""
    communitiesResult = "Communities (Undirected):\n"
    graphForCommunities = nx.Graph()
    bestPartition = dict()
    listPartition = dict()
    networkModularity = 0
    everyCommunity = dict()
    listCommunities = []
    colors = []
    cmap = []
    between = []
    eachNodes = []
    eachColor = []
    positions = []
    name = CleanNetworkFilename(fileName)
    counter = 0
    counter2 = 0
    # get undir graph integer-weighted
    graphForCommunities = ParseFileToNetwork(fileName, "-u", "community")
    # get communities
    bestPartition = community.best_partition(graphForCommunities)
    networkModularity = community.modularity(bestPartition, graphForCommunities)
    listPartition = list(bestPartition.keys())
    for counter in range(len(listPartition)):
        if(not(bestPartition[listPartition[counter]] in everyCommunity)):
            everyCommunity[bestPartition[listPartition[counter]]] = []
            everyCommunity[bestPartition[listPartition[counter]]].append(listPartition[counter])
        else:
            everyCommunity[bestPartition[listPartition[counter]]].append(listPartition[counter])
    communitiesResult = communitiesResult + "- Modularity: " + str(networkModularity) + "\n"
    listCommunities = list(everyCommunity.keys())
    communitiesResult = communitiesResult + "- Number of Communities: " + str(len(listCommunities)) + "\n"
    for counter in range(len(listCommunities)):
        communitiesResult=communitiesResult+"- Number of nodes in community "+str(counter+1)+": "+str(len(everyCommunity[listCommunities[counter]]))+"\n"
        communitiesResult = communitiesResult + "- Nodes in C_" + str(counter + 1) + ": " + ",".join(everyCommunity[listCommunities[counter]]) + "\n"
    # plot message
    print("\t- Drawing networks ...")
    # plot graph with communities
    cmap = plt.cm.get_cmap('gist_rainbow')
    between = arange(0.01, 1, float(0.99/float(len(listCommunities))))
    for counter in range(len(listCommunities)):
        for counter2 in range(len(everyCommunity[listCommunities[counter]])):
            eachNodes.append(everyCommunity[listCommunities[counter]][counter2])
            eachColor.append(cmap(between[counter]))
    positions = nx.spring_layout(graphForCommunities)
    nx.draw_networkx(graphForCommunities, with_labels = False, pos = positions, nodelist = eachNodes, node_color = eachColor, node_size = 20, width = 0.4)
    plt.axis("off")
    plt.title(name + " communities")
    plt.tight_layout()
    plt.savefig(  name + "_community_network.png" , dpi=300)
    
    plt.close()
    if(someNetwork.is_directed()):
        nx.draw_networkx(someNetwork, with_labels = False, pos = positions, node_size = 20, width = 0.4)
        plt.axis("off")
        plt.title(name)
        plt.savefig(  name + "_directed_network.png" , dpi=300) 
        plt.tiight_layout()
        plt.close()
    else:
        nx.draw_networkx(someNetwork, with_labels = False, pos = positions, node_size = 20, width = 0.4)
        plt.axis("off")
        plt.title(name)
        plt.savefig( name + "_undirected_network.png" , dpi=300) 
        plt.tight_layout()
        plt.close()
    # end of function
    communitiesFinalResult = "\n\n" + communitiesResult + "\n\n"
    return(communitiesFinalResult)

# Function: analize degree distribution ----------------------------------------------------------------------------------------------------------------------
def AnalizeDegreeDistribution(someNetwork, typeNW, fileName):
    # function message
    print("\t- Obtaining degree distribution ...")
    # variables
    DegDistributionFinalResult = ""
    DegDistributionResult = "Degree Distribution:\n"
    name = CleanNetworkFilename(fileName)
    degMapping = dict()
    inDegMapping = dict()
    outDegMapping = dict()
    maxDegree = 0
    inMaxDegree = 0
    outMaxDegree = 0
    minDegree = 0
    inMinDegree = 0
    outMinDegree = 0
    degVector = dict()
    inDegVector = dict()
    outDegVector = dict()
    distVector = []
    inDistVector = []
    outDistVector = []
    bars = []
    counter = 0
    strDistVector = []
    strInDistVector = []
    strOutDistVector = []
    # get degree distribution for undirected
    if(typeNW == "-u"):
        degMapping = dict(someNetwork.degree())
        degVector = list(degMapping.values())
        degVector.sort()
        minDegree = degVector[0]
        degVector.reverse()
        maxDegree = degVector[0]
        # drawing message
        print("\t- Drawing Degree Distribution ...")
        # draw degree distribution
        for counter in range(maxDegree + 1):
            distVector.append(0)
        for counter in range(len(degVector)):
            distVector[degVector[counter]] = distVector[degVector[counter]] + 1
        bars = range(len(distVector))
        plt.bar(bars, distVector, bottom = 0, width = 0.99, align = "center", color = "m")
        plt.xticks(range(0, maxDegree + 1, int(math.ceil(maxDegree*0.05))))
        plt.xlabel("Degree")
        plt.ylabel("Nodes with each degree")
        plt.title(name + " degree distribution")
        plt.tight_layout()
        plt.savefig(  name + "_distribution_degree.png",dpi=300)
        plt.close()
        for counter in range(len(distVector)):
            strDistVector.append(str(distVector[counter]))
        DegDistributionResult = DegDistributionResult + "- Degree distribution:" + ",".join(list(strDistVector))
    # get degree distribution for directed
    if(typeNW == "-d"):
        inDegMapping = dict(someNetwork.in_degree())
        outDegMapping = dict(someNetwork.out_degree())
        inDegVector = list(inDegMapping.values())
        outDegVector = list(outDegMapping.values())
        inDegVector.sort()
        outDegVector.sort()
        inMinDegree = inDegVector[0]
        outMinDegree = outDegVector[0]
        inDegVector.reverse()
        outDegVector.reverse()
        inMaxDegree = inDegVector[0]
        outMaxDegree = outDegVector[0]
        # drawing message
        print("\t- Drawing in degree distribution ...")
        # draw in degree distribution
        for counter in range(inMaxDegree + 1):
            inDistVector.append(0)
        for counter in range(len(inDegVector)):
            inDistVector[inDegVector[counter]] = inDistVector[inDegVector[counter]] + 1
        bars = range(len(inDistVector))
        plt.bar(bars, inDistVector, bottom = 0, width = 0.99, align = "center", color = "b")
        plt.xticks(range(0, inMaxDegree + 1, int(math.ceil(inMaxDegree*0.05))))
        plt.xlabel("In Degree")
        plt.ylabel("Nodes with each in-degree")
        plt.title(name + " in-degree distribution")
        plt.tight_layout()
        plt.savefig( name + "distribution_indegree.png", dpi=300)
       
        plt.close()
        for counter in range(len(inDistVector)):
            strInDistVector.append(str(inDistVector[counter]))
        DegDistributionResult = DegDistributionResult + "- In Degree distribution:" + ",".join(list(strInDistVector)) + "\n"
        # drawing message
        print("\t- Drawing out degree distribution ...")
        # draw out degree distribution
        for counter in range(outMaxDegree + 1):
            outDistVector.append(0)
        for counter in range(len(outDegVector)):
            outDistVector[outDegVector[counter]] = outDistVector[outDegVector[counter]] + 1
        bars = range(len(outDistVector))
        plt.bar(bars, outDistVector, bottom = 0, width = 0.99, align = "center", color = "r")
        plt.xticks(range(0, outMaxDegree + 1, int(math.ceil(outMaxDegree*0.05))))
        plt.xlabel("Out Degree")
        plt.ylabel("Nodes with each out-degree")
        plt.title(name + " out-degree distribution")
        plt.tight_layout()
        plt.savefig( name + "distribution_outdegree.png", dpi=300)
        plt.close()
        for counter in range(len(outDistVector)):
            strOutDistVector.append(str(outDistVector[counter]))
        DegDistributionResult = DegDistributionResult + "- Out Degree distribution:" + ",".join(list(strOutDistVector))
    # end of function
    DegDistributionFinalResult = "\n\n" + DegDistributionResult + "\n\n"
    return(DegDistributionFinalResult)

# Function: analize network -----------------------------------------------
def AnalizeNetwork(analysisNetwork, networkFileName, typeNW):
    # variables
    totalResults = ""
    totalResultsFile = None
    name = CleanNetworkFilename(networkFileName)
    # start analysis
    totalResults = totalResults + AnalizeOrder(analysisNetwork)
    totalResults = totalResults + AnalizeSize(analysisNetwork)
    totalResults = totalResults + AnalizeDiameter(analysisNetwork)
    totalResults = totalResults + AnalizeRadius(analysisNetwork)
    totalResults = totalResults + AnalizeDensity(analysisNetwork)
    totalResults = totalResults + AnalizeMeanDegree(analysisNetwork, typeNW)
    totalResults = totalResults + AnalizeMaxDegree(analysisNetwork)
    totalResults = totalResults + AnalizeClustCoeff(analysisNetwork, typeNW)
    totalResults = totalResults + AnalizeMaxHubs(analysisNetwork, typeNW)
    totalResults = totalResults + AnalizeCutVertices(analysisNetwork)
    totalResults = totalResults + AnalizeCutEdges(analysisNetwork)
    totalResults = totalResults + AnalizeConnComps(analysisNetwork)
    totalResults = totalResults + AnalizeStronConnComps(analysisNetwork, typeNW)
    totalResults = totalResults + AnalizeMaxCliques(analysisNetwork)
    totalResults = totalResults + AnalizeCycleBasis(analysisNetwork)
    totalResults = totalResults + AnalizeIndependentSets(analysisNetwork, 5)
    totalResults = totalResults + AnalizeDegreeDistribution(analysisNetwork, typeNW, networkFileName)
    totalResults = totalResults + AnalizeCommunitiesAndDrawings(networkFileName, analysisNetwork)
    # open results file
    if(typeNW == "-u"):
        totalResultsFile = open( name,"_undirected.txt" "w")
    if(typeNW == "-d"):
        totalResultsFile = open( name, "_directed.txt " "w")
    # print results
    totalResultsFile.write(totalResults)
    totalResultsFile.close()
    # end of function

# Function: analize random network -------------------------------------
def RandomNetworkAnalysis(someNetwork):
    # variables
    state = "UNDER CONSTRUCTION"

# Function: check for correct input -----------------------------------
def CheckInputStr(inList):
    # variables
    totalIn = len(inList)
    right = 1
    # check argv right format
    if(totalIn%2==1):
        for i in range(1, totalIn):
            if(i%2 == 0):
                if(".txt" in inList[i]):
                    continue
                else:
                    right = 0
                    break
            if(i%2 == 1):
                if((inList[i]=="-u") or (inList[i]=="-d")):
                    continue
                else:
                    right = 0
                    break
    else:
        right = 0
    # end of function
    return(right)

# MAIN ##############################################################

# check if command line input is given correctly
if(CheckInputStr(argv) == 1):
    # Analize every network in command line
    for networkCounter in range(1, len(argv), 2):
        networkType = argv[networkCounter]
        if(networkType == "-u"):
            print("\n> Analysing Network:\t" + argv[networkCounter] + "\ttype:\tUndirected")
            originalUndirNetwork = ParseFileToNetwork(argv[networkCounter + 1], networkType, "else")
            AnalizeNetwork(originalUndirNetwork, argv[networkCounter + 1], networkType)
            print("\t- Finished analysis for this network")
        if(networkType == "-d"):
            print("\n> Analysing Network:\t" + argv[networkCounter] + "\ttype:\tDirected")
            originalDirecNetwork = ParseFileToNetwork(argv[networkCounter + 1], networkType, "else")
            AnalizeNetwork(originalDirecNetwork, argv[networkCounter + 1], networkType)
            print("\t- Finished analysis for this network")
    print("\n\n--- Finished Analyzing Every Given Network ---\n\n")
else:
    print("\n\n\t :O Sorry, Wrong Input.\n\n")

# SUPLEMENTARY ####################################################

# Core function code -----------------------------------------------------------------------------------------------------------------------------------------
#def AnalizeInvariant(someNetwork):
#    # function message
#    print("\t- Analysing Invarian ...")
#    # variables
#    invariantFinalResult = ""
#    invariantResult = "Invariant:\n"
#    # get invariant
#    # end of function
#    invariantFinalResult = "\n\n" + sizeResult + "\n\n//////////////////////////////////////////////////////////////////////"
#    return(invariantFinalResult)

#########################################################################
