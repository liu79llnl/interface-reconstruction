import numpy as np
import math
import vtk
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from geoms import getDistance, getArea, mergePolys, getPolyIntersectArea, getPolyLineArea, getPolyLineIntersects, lineIntersect, getCentroid
from linear_facet import getLinearFacet
from circular_facet import getCircleIntersectArea, getCircleCircleIntersects, getArcFacet, getArcFacetNewton, getCircleLineIntersects2, getCenter
from corner_facet import getPolyCornerArea, getPolyCurvedCornerArea, getCurvedCornerFacet
from mesh import makeCartesianGrid, makeQuadGrid, makeConcaveGrid, makeFineCartesianGrid
from vtkplots import plotQuadGrid, plotFacets

#opolys = grid of quads
#areas = grid of area fractions per quad
def merge(opolys, areas):
    #Merging
    dirs = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    mergedpolys = [] #coordinates of polys to be merged
    predmergedpolys = [[[] for _ in range(len(opolys[0]))] for _ in range(len(opolys))]
    neighborsmergedpolys = [[[] for _ in range(len(opolys[0]))] for _ in range(len(opolys))]

    #Compute initial neighbors
    for x in range(len(opolys)):
        for y in range(len(opolys[0])):
            if areas[x][y] == 0 or areas[x][y] == 1:
                predmergedpolys[x][y] = None
                neighborsmergedpolys[x][y] = None
                continue
            #Compute number of neighbors
            gridsquare = [x, y]
            predmergedpolys[x][y] = [[x, y]]
            for direction in dirs:
                #Check directions in counterclockwise order
                testsquare = [gridsquare[0] + direction[0], gridsquare[1] + direction[1]]
                if testsquare[0] < len(opolys) and testsquare[0] >= 0 and testsquare[1] < len(opolys[0]) and testsquare[1] >= 0 and areas[testsquare[0]][testsquare[1]] < 1 and areas[testsquare[0]][testsquare[1]] > 0:
                    #This direction is a partial area neighbor
                    neighborsmergedpolys[x][y].append(testsquare)
    
    #By here, predmergedpolys[x][y] = None if full, [[x, y]] otherwise. neighborsmergedpolys[x][y] = None if full, [neighbors] otherwise.
    #Merge squares with 1 or 3+ neighbors
    for x in range(len(opolys)):
        for y in range(len(opolys[0])):
            gridneighbors = neighborsmergedpolys[x][y]
            if gridneighbors is None:
                #Full cell, skip
                continue
            #Casework on number of neighbors
            assert len(gridneighbors) > 0, "Mixed cell with no mixed neighbors: increase fineness of resolution!"
            if len(gridneighbors) == 1 and not(x == 0 or y == 0 or x == len(opolys)-1 or y == len(opolys[0])-1):
                #Long thin strip, probably corner here
                mergesquares = [[x, y]]
                curx = gridneighbors[0][0]
                cury = gridneighbors[0][1]
                prevx = x
                prevy = y
                continueSpike = True
                while continueSpike:
                    mergesquares.append(opolys[curx][cury])
                    if len(neighborsmergedpolys[curx][cury]) >= 3:
                        continueSpike = False
                        #Update merged cells
                        mergeneighbors = []
                        for mergeneighbor in neighborsmergedpolys[curx][cury]:
                            if mergeneighbor != [prevx, prevy]:
                                mergeneighbors.append(mergeneighbor)
                        mergecells = []
                        for mergesquare in mergesquares:
                            for mergecell in predmergedpolys[mergesquare[0]][mergesquare[1]]:
                                if mergecell not in mergecells:
                                    mergecells.append(mergecell)
                        for mergecell in mergecells:
                            predmergedpolys[mergecell[0]][mergecell[1]] = mergecells
                            neighborsmergedpolys[mergecell[0]][mergecell[1]] = mergeneighbors
                    else:
                        assert len(neighborsmergedpolys[curx][cury]) == 2, "Mixed cells form one cell wide strip: increase fineness of resolution!"
                        if neighborsmergedpolys[curx][cury][0] == [prevx, prevy]:
                            prevx = curx
                            prevy = cury
                            cursquare = neighborsmergedpolys[curx][cury][1]
                            curx = cursquare[0]
                            cury = cursquare[1]
                        else:
                            prevx = curx
                            prevy = cury
                            cursquare = neighborsmergedpolys[curx][cury][0]
                            curx = cursquare[0]
                            cury = cursquare[1]
            else:
                #3+ neighbors
                #invariant: mergeneighbors contains cells to be merged, whose neighbors may not have been explored yet; gridneighbors are the cells that need to be explored
                mergecells = []
                for mergecell in predmergedpolys[x][y]:
                    mergecells.append(mergecell)
                mergeneighbors = []
                for gridneighbor in gridneighbors:
                    for gridcell in predmergedpolys[gridneighbor[0]][gridneighbor[1]]:
                        mergecells.append(gridcell)
                    mergeneighbors.append(gridneighbor)
                while len(mergeneighbors) > 0:
                    testsquare = mergeneighbors[0]
                    mergeneighbors.pop(0)
                    testneighbors = neighborsmergedpolys[testsquare[0]][testsquare[1]]
                    if len(testneighbors) >= 3:
                        for testneighbor in testneighbors:
                            if testneighbor not in mergecells:
                                for mergecell in predmergedpolys[testneighbor[0]][testneighbor[1]]:
                                    mergecells.append(mergecell)
                                mergeneighbors.append(testneighbor)
                    elif len(testneighbors) == 2:
                        for testneighbor in testneighbors:
                            if testneighbor not in mergecells:
                                dircount = 0
                                for test2neighbor in neighborsmergedpolys[testneighbor[0]][testneighbor[1]]:
                                    if test2neighbor in mergecells:
                                        dircount += 1
                                if dircount >= 2:
                                    for mergecell in predmergedpolys[testneighbor[0]][testneighbor[1]]:
                                        mergecells.append(mergecell)
                                    mergeneighbors.append(testneighbor)
                #mergecells consists of cells to be merged + the two neighbors (not to be merged)
                #TODO: current code greedily merges all ambiguous cells until ambiguities gone. This can create large, narrow clusters of merged cells, which ideally could be broken into smaller clusters.
                mergeneighbors = []
                for mergecell in mergecells:
                    nummergeneighbors = 0
                    numallneighbors = len(neighborsmergedpolys[mergecell[0]][mergecell[1]])
                    for testneighbor in neighborsmergedpolys[mergecell[0]][mergecell[1]]:
                        if testneighbor in mergecells:
                            nummergeneighbors += 1
                    if nummergeneighbors < 2 and numallneighbors > 1:
                        #Endpoint
                        mergeneighbors.append(mergecell)
                for mergeneighbor in mergeneighbors:
                    mergecells.remove(mergeneighbor)
                #Update merged cells
                for mergecell in mergecells:
                    predmergedpolys[mergecell[0]][mergecell[1]] = mergecells
                    neighborsmergedpolys[mergecell[0]][mergecell[1]] = mergeneighbors
    
    mergedpolyindices = [[None for _ in range(len(opolys[0]))] for _ in range(len(opolys))]
    mergedpolyinfos = []
    mergedcoords = [] #the coordinates of the fully merged poly
    mergedareafractions = [] #area fractions of fully merged poly

    for x in range(len(opolys)):
        for y in range(len(opolys[0])):
            if len(predmergedpolys[x][y]) == 0:
                continue
            elif len(predmergedpolys[x][y]) == 1:
                #x, y
                mergedpolyindices[x][y] = len(mergedpolyinfos)
                mergedpolyinfos.append([predmergedpolys[x][y], neighborsmergedpolys[x][y]])
                mergedcoords.append(opolys[x][y])
                mergedareafractions.append(areas[x][y])
            else:
                #cluster
                mergedpolyinfo = [predmergedpolys[x][y], neighborsmergedpolys[x][y]]
                if mergedpolyinfo not in mergedpolyinfos:
                    mergecells = predmergedpolys[x][y]
                    for mergecell in mergecells:
                        mergedpolyindices[mergecell[0]][mergecell[1]] = len(mergedpolyinfos)
                    mergedpolyinfos.append(mergedpolyinfo)
                    boundary = mergecells[0]
                    for mergecelli in range(1, len(mergecells)):
                        boundary = mergePolys(boundary, opolys[mergecells[mergecelli][0]][mergecells[mergecelli][1]])
                    mergedcoords.append(boundary)
                    boundaryarea = sum(list(map(lambda x : areas[x[0]][x[1]]*getArea(opolys[x[0]][x[1]]), mergecells)))/sum(list(map(lambda x : getArea(opolys[x[0]][x[1]]), mergecells)))
                    mergedareafractions.append(boundaryarea)

    return mergedpolyindices, mergedpolyinfos, mergedcoords, mergedareafractions