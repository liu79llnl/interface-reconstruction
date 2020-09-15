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
            if gridneighbors is None or len(gridneighbors) == 2:
                #Full cell or unambiguous cell, skip
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
            if predmergedpolys[x][y] is None or len(predmergedpolys[x][y]) == 0:
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
                    boundary = opolys[mergecells[0][0]][mergecells[0][1]]
                    for mergecelli in range(1, len(mergecells)):
                        boundary = mergePolys(boundary, opolys[mergecells[mergecelli][0]][mergecells[mergecelli][1]])
                    mergedcoords.append(boundary)
                    boundaryarea = sum(list(map(lambda x : areas[x[0]][x[1]]*getArea(opolys[x[0]][x[1]]), mergecells)))/sum(list(map(lambda x : getArea(opolys[x[0]][x[1]]), mergecells)))
                    mergedareafractions.append(boundaryarea)

    return mergedpolyindices, mergedpolyinfos, mergedcoords, mergedareafractions

#mergedpolyindices: xy-grid of merged polygon indices: if mergedpolyindices[x][y] = i, mergedpolyinfos[i] is the merged polygon at [x, y]
#mergedpolyinfos: list of merged polygons: [[list of [x, y]s to be merged], [[neighbor1], [neighbor2]]]
#mergedcoords: list of coords for merged polygons in same order as mergedpolyinfos
#mergedareafractions: list of area fractions for merged polygons in same order as mergedpolyinfos

#mergedpolyindices, mergedpolyinfos, mergedcoords, mergedareafractions: from merging function
#opolys, areas: original volume fractions, used to refine linear facets
#threshold: tolerance for errors in volume fractions
def makeFacets(mergedpolyindices, mergedpolyinfos, mergedcoords, mergedareafractions, opolys, areas, threshold, makeGapless):
    #facets
    predfacets = [[None] * len(opolys[0]) for _ in range(len(opolys))]
    facetsunique = []

    #orientations
    predorientations = [[None] * len(opolys[0]) for _ in range(len(opolys))]
    mergedorientations = [None] * len(mergedpolyinfos)
    
    #compute facet orientations
    for y in range(len(opolys[0])):
        for x in range(len(opolys)):
            if mergedpolyindices[x][y] is not None and mergedorientations[mergedpolyindices[x][y]] == None:
                #finding path, starting with bottom left most point
                #orientation of path: outer boundary = True, inner boundary = False
                if (y > 0 and areas[x][y-1] == 1) or (x > 0 and areas[x-1][y] == 1):
                    curOrientation = False
                else:
                    curOrientation = True
                
                path = []
                curmergedpolyindex = mergedpolyindices[x][y]
                curmergedpoly = mergedpolyinfos[curmergedpolyindex]
                if curmergedpoly[1][0][0] < curmergedpoly[1][1][0]:
                    nextmergedpolyindex = mergedpolyindices[curmergedpoly[1][1][0]][curmergedpoly[1][1][1]]
                elif curmergedpoly[1][0][0] > curmergedpoly[1][1][0]:
                    nextmergedpolyindex = mergedpolyindices[curmergedpoly[1][0][0]][curmergedpoly[1][0][1]]
                else:
                    if curmergedpoly[1][0][1] < curmergedpoly[1][1][1]:
                        nextmergedpolyindex = mergedpolyindices[curmergedpoly[1][0][0]][curmergedpoly[1][0][1]]
                    else:
                        nextmergedpolyindex = mergedpolyindices[curmergedpoly[1][1][0]][curmergedpoly[1][1][1]]
                path.append(curmergedpolyindex)
                mergedorientations[curmergedpolyindex] = curOrientation
                for point in curmergedpoly[0]:
                    predorientations[point[0]][point[1]] = curOrientation
                prevmergedpolyindex = curmergedpolyindex
                while nextmergedpolyindex != path[0]:
                    #Here: handle single neighbor cases
                    if nextmergedpolyindex is not None:
                        curmergedpolyindex = nextmergedpolyindex
                    else:
                        print(mergedpolyinfos[curmergedpolyindex])
                        print(1/0)
                    curmergedpoly = mergedpolyinfos[curmergedpolyindex]
                    if curmergedpoly[1][0] in mergedpolyinfos[prevmergedpolyindex][0]:
                        nextmergedpolyindex = mergedpolyindices[curmergedpoly[1][1][0]][curmergedpoly[1][1][1]]
                    else:
                        assert curmergedpoly[1][1] in mergedpolyinfos[prevmergedpolyindex][0], print(mergedpolyinfos[prevmergedpolyindex])
                        nextmergedpolyindex = mergedpolyindices[curmergedpoly[1][0][0]][curmergedpoly[1][0][1]]
                    path.append(curmergedpolyindex)
                    mergedorientations[curmergedpolyindex] = curOrientation
                    for point in curmergedpoly[0]:
                        predorientations[point[0]][point[1]] = curOrientation
                    prevmergedpolyindex = curmergedpolyindex
                    
                assert len(path) >= 2
                path.append(path[0])
                path.append(path[1])
                
                facetfitted = [None] * (len(path)-2)
                
                #Make linear facets
                for pathelement in range(1, len(path)-1):
                    #previ, curi, nexti are indices of mergedpolyinfos
                    curi = path[pathelement]
                    #invert path if needed
                    if mergedorientations[curi]:
                        previ = path[pathelement-1]
                        nexti = path[pathelement+1]
                    else:
                        previ = path[pathelement+1]
                        nexti = path[pathelement-1]
                    prevpolygon = mergedcoords[previ]
                    prevpolygonarea = mergedareafractions[previ]
                    nextpolygon = mergedcoords[nexti]
                    nextpolygonarea = mergedareafractions[nexti]
                    try:
                        facetline1, facetline2 = getLinearFacet(prevpolygon, nextpolygon, prevpolygonarea, nextpolygonarea, threshold)
                        if abs(mergedareafractions[curi] - getPolyLineArea(mergedcoords[curi], facetline1, facetline2)/getArea(mergedcoords[curi])) > threshold:
                            #maybe excess merging reduced resolution: try immediate neighbors
                            if mergedpolyinfos[curi][1][0] in mergedpolyinfos[previ][0]:
                                prevneighbor = mergedpolyinfos[curi][1][0]
                                nextneighbor = mergedpolyinfos[curi][1][1]
                            else:
                                assert mergedpolyinfos[curi][1][1] in mergedpolyinfos[previ][0]
                                prevneighbor = mergedpolyinfos[curi][1][1]
                                nextneighbor = mergedpolyinfos[curi][1][0]
                            prevpolygon = opolys[prevneighbor[0]][prevneighbor[1]]
                            nextpolygon = opolys[nextneighbor[0]][nextneighbor[1]]
                            prevpolygonarea = areas[prevneighbor[0]][prevneighbor[1]]
                            nextpolygonarea = areas[nextneighbor[0]][nextneighbor[1]]
                            facetline1, facetline2 = getLinearFacet(prevpolygon, nextpolygon, prevpolygonarea, nextpolygonarea, threshold)
                        if abs(mergedareafractions[curi] - getPolyLineArea(mergedcoords[curi], facetline1, facetline2)/getArea(mergedcoords[curi])) < linearerrorthreshold:
                            intersects = getPolyLineIntersects(mergedcoords[curi], facetline1, facetline2)
                            if len(intersects) > 0:
                                facetfitted[pathelement-1] = ['linear', intersects]
                        else:
                            continue
                    except:
                        print("Failed linear facet")
                            
                #Make corners
                #Find closest linear facets to left/right of corner poly
                lefts = [None] * (len(facetfitted))
                rights = [None] * (len(facetfitted))
                curlinearfacet = None
                for pathelement in range(1, 2*len(facetfitted)-1):
                    if facetfitted[(pathelement) % len(facetfitted)] is None:
                        lefts[(pathelement) % len(facetfitted)] = curlinearfacet
                    else:
                        curlinearfacet = facetfitted[(pathelement) % len(facetfitted)][1]
                for pathelement in range(1, 2*len(facetfitted)-1):
                    if facetfitted[(len(facetfitted)-1-pathelement) % len(facetfitted)] is None:
                        rights[(len(facetfitted)-1-pathelement) % len(facetfitted)] = curlinearfacet
                    else:
                        curlinearfacet = facetfitted[(len(facetfitted)-1-pathelement) % len(facetfitted)][1]

                #Try to make corners
                #here pathelement is 1 less than pathelement from linear facet fitting
                for pathelement in range(len(facetfitted)):
                    if facetfitted[pathelement] is None and lefts[pathelement] is not None and rights[pathelement] is not None and lefts[pathelement][0] != rights[pathelement][0] and lefts[pathelement][1] != rights[pathelement][1]:
                        #invert path if needed
                        if mergedorientations[path[pathelement]]:
                            prevFacet = lefts[pathelement]
                            nextFacet = rights[pathelement]
                        else:
                            prevFacet = rights[pathelement]
                            nextFacet = lefts[pathelement]
                        corner = lineIntersect(prevFacet[0], prevFacet[1], nextFacet[0], nextFacet[1])
                        corner = [prevFacet[0], corner, nextFacet[1]]
                        cornerareafraction = getPolyCornerArea(mergedcoords[path[pathelement+1]], corner[0], corner[1], corner[2])/getArea(mergedcoords[path[pathelement+1]])
                        if abs(cornerareafraction - mergedareafractions[path[pathelement+1]]) < cornererrorthreshold:
                            facetfitted[pathelement] = ['corner', corner]
                        else:
                            continue

                #Try to make arc facets for the remaining ones
                for pathelement in range(1, len(path)-1):
                    if facetfitted[pathelement-1] is None:
                        #previ, curi, nexti are indices of mergedpolyinfos
                        curi = path[pathelement]
                        #invert path if needed
                        if mergedorientations[grid2]:
                            previ = path[pathelement-1]
                            nexti = path[pathelement+1]
                        else:
                            previ = path[pathelement+1]
                            nexti = path[pathelement-1]
                        prevpolygon = mergedcoords[previ]
                        prevpolygonarea = mergedareafractions[previ]
                        curpolygon = mergedcoords[curi]
                        curpolygonarea = mergedareafractions[curi]
                        nextpolygon = mergedcoords[nexti]
                        nextpolygonarea = mergedareafractions[nexti]
                        try:
                            arccenter, arcradius, arcintersects = getArcFacet(prevpolygon, curpolygon, nextpolygon, prevpolygonarea, curpolygonarea, nextpolygonarea, threshold)
                            if len(arcintersects) == 0:
                                print("getArcFacet({}, {}, {}, {}, {}, {}, {})".format(prevpolygon, curpolygon, nextpolygon, prevpolygonarea, curpolygonarea, nextpolygonarea, threshold))
                            else:
                                facetfitted[pathelement-1] = ['arc', arccenter, arcradius, arcintersects]
                        except:
                            print("getArcFacet({}, {}, {}, {}, {}, {}, {})".format(prevpolygon, curpolygon, nextpolygon, prevpolygonarea, curpolygonarea, nextpolygonarea, threshold))
                            continue

                #Retry failed circle facets with circle facet adjacent guess
                for i in range(len(facetfitted)):
                    if facetfitted[i] is None or (facetfitted[i][0] == 'arc' and len(facetfitted[i][3]) == 0):
                        #previ, curi, nexti are indices of mergedpolyinfos
                        curi = path[i+1]
                        #invert path if needed
                        if mergedorientations[curi]:
                            previ = path[i]
                            nexti = path[i+2]
                        else:
                            previ = path[i+2]
                            nexti = path[i]
                        prevpolygon = mergedcoords[previ]
                        prevpolygonarea = mergedareafractions[previ]
                        curpolygon = mergedcoords[curi]
                        curpolygonarea = mergedareafractions[curi]
                        nextpolygon = mergedcoords[nexti]
                        nextpolygonarea = mergedareafractions[nexti]
                        #Try previous to see if circle facet
                        if facetfitted[(i-1) % len(facetfitted)] is not None and facetfitted[(i-1) % len(facetfitted)][0] == 'arc':
                            gcenterx = facetfitted[(i-1) % len(facetfitted)][1][0]
                            gcentery = facetfitted[(i-1) % len(facetfitted)][1][1]
                            gradius = facetfitted[(i-1) % len(facetfitted)][2]
                        #Try next to see if circle facet
                        elif facetfitted[(i+1) % len(facetfitted)] is not None and facetfitted[(i+1) % len(facetfitted)][0] == 'arc':
                            gcenterx = facetfitted[(i+1) % len(facetfitted)][1][0]
                            gcentery = facetfitted[(i+1) % len(facetfitted)][1][1]
                            gradius = facetfitted[(i+1) % len(facetfitted)][2]
                        try:
                            arccenter, arcradius, arcintersects = getArcFacetNewton(prevpolygon, curpolygon, nextpolygon, prevpolygonarea, curpolygonarea, nextpolygonarea, threshold, gcenterx, gcentery, gradius)
                            facetfitted[i] = ['arc', arccenter, arcradius, arcintersects]
                        except:
                            #Circle facet algorithm with Newton's method failed or no nearby circular facets
                            print("getArcFacetNewton({}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(prevpolygon, curpolygon, nextpolygon, prevpolygonarea, curpolygonarea, nextpolygonarea, threshold, gcenterx, gcentery, gradius))
                            #Try to fit a curved corner
                            #Find leftward arc neighbor
                            prevarc = (i-1) % len(facetfitted)
                            #Go until you find two arc neighboring arc facets that have similar curvature
                            while prevarc != i and (facetfitted[prevarc] is None or facetfitted[(prevarc-1) % len(facetfitted)] is None or (facetfitted[prevarc][0] != 'arc' and facetfitted[prevarc][0] != 'linear') or (facetfitted[prevarc][0] == 'arc' and facetfitted[(prevarc-1) % len(facetfitted)][0] == 'arc' and abs(facetfitted[prevarc][2]-facetfitted[(prevarc-1) % len(facetfitted)][2]) > similararcthreshold)):
                                prevarc = (prevarc-1) % len(facetfitted)
                            if prevarc == i:
                                #No neighbors
                                continue
                            nextarc = prevarc+1
                            #Go until you find two arc neighboring arc facets that have similar curvature
                            while nextarc != prevarc+len(facetfitted) and (facetfitted[nextarc % len(facetfitted)] is None or facetfitted[(nextarc+1) % len(facetfitted)] is None or (facetfitted[nextarc % len(facetfitted)][0] != 'arc' and facetfitted[nextarc % len(facetfitted)][0] != 'linear') or (facetfitted[nextarc % len(facetfitted)][0] == 'arc' and facetfitted[(nextarc+1) % len(facetfitted)][0] == 'arc' and abs(facetfitted[nextarc % len(facetfitted)][2]-facetfitted[(nextarc+1) % len(facetfitted)][2]) > similararcthreshold)):
                                nextarc += 1
                            if nextarc % len(facetfitted) == prevarc:
                                #Only one neighbor
                                continue
                            nextarc = nextarc % len(facetfitted)
                            prevcornerfacetside = facetfitted[prevarc]
                            nextcornerfacetside = facetfitted[nextarc]
                            if prevcornerfacetside[0] == 'linear' and nextcornerfacetside[0] == 'arc':
                                prevcenter = None
                                prevradius = None
                                curvedcorner1 = prevcornerfacetside[1][0]
                                nextcenter = nextcornerfacetside[1]
                                nextradius = nextcornerfacetside[2]
                                curvedcorner2 = nextcornerfacetside[3][0]

                                curvedcornerpoint = getCircleLineIntersects2(prevcornerfacetside[1][0], prevcornerfacetside[1][1], nextcenter, nextradius)
                                if len(curvedcornerpoint) == 1:
                                    curvedcornerpoint = curvedcornerpoint[0]
                                elif len(curvedcornerpoint) > 1:
                                    if getDistance(curpolygon[0], curvedcornerpoint[0]) >= getDistance(curpolygon[0], curvedcornerpoint[1]):
                                        curvedcornerpoint = curvedcornerpoint[1]
                                    else:
                                        curvedcornerpoint = curvedcornerpoint[0]
                            elif prevcornerfacetside[0] == 'arc' and nextcornerfacetside[0] == 'linear':
                                prevcenter = prevcornerfacetside[1]
                                prevradius = prevcornerfacetside[2]
                                curvedcorner1 = prevcornerfacetside[3][0]
                                nextcenter = None
                                nextradius = None
                                curvedcorner2 = nextcornerfacetside[1][1]

                                curvedcornerpoint = getCircleLineIntersects2(nextcornerfacetside[1][0], nextcornerfacetside[1][1], prevcenter, prevradius)
                                if len(curvedcornerpoint) == 1:
                                    curvedcornerpoint = curvedcornerpoint[0]
                                elif len(curvedcornerpoint) > 1:
                                    if getDistance(curpolygon[0], curvedcornerpoint[0]) >= getDistance(curpolygon[0], curvedcornerpoint[1]):
                                        curvedcornerpoint = curvedcornerpoint[1]
                                    else:
                                        curvedcornerpoint = curvedcornerpoint[0]
                            elif prevcornerfacetside[0] == 'arc' and nextcornerfacetside[0] == 'arc':
                                prevcenter = prevcornerfacetside[1]
                                prevradius = prevcornerfacetside[2]
                                curvedcorner1 = prevcornerfacetside[3][-1]
                                nextcenter = nextcornerfacetside[1]
                                nextradius = nextcornerfacetside[2]
                                curvedcorner2 = nextcornerfacetside[3][0]

                                curvedcornerpoint = getCircleCircleIntersects(prevcenter, nextcenter, prevradius, nextradius)
                                if getDistance(curpolygon[0], curvedcornerpoint[0]) >= getDistance(curpolygon[0], curvedcornerpoint[1]):
                                    curvedcornerpoint = curvedcornerpoint[1]
                                else:
                                    curvedcornerpoint = curvedcornerpoint[0]
                            else:
                                #two linears
                                prevradius = None
                                curvedcorner1 = prevcornerfacetside[1][0]
                                nextradius = None
                                curvedcorner2 = nextcornerfacetside[1][1]
                                curvedcornerpoint = lineIntersect(prevcornerfacetside[1][0], prevcornerfacetside[1][1], nextcornerfacetside[1][0], nextcornerfacetside[1][1])
                            
                            #Use all polys in between two nearest successful facets
                            curvedcornergridindices = list(map(lambda x : path[(x % len(facetfitted))+1], range(prevarc+1, nextarc)))
                            curvedcornerpolys = list(map(lambda x : mergedcoords[x], curvedcornergridindices))
                            curvedcornerareas = list(map(lambda x : mergedareafractions[x], curvedcornergridindices))
                            try:
                                if mergedorientations[path[i+1]]:
                                    curvedcorner = getCurvedCornerFacet(curvedcornerpolys, curvedcornerareas, curvedcorner1, curvedcornerpoint, curvedcorner2, prevradius, nextradius, threshold, curvedconvthreshold, curvednewtonfactor)
                                else:
                                    curvedcorner = getCurvedCornerFacet(curvedcornerpolys, curvedcornerareas, curvedcorner2, curvedcornerpoint, curvedcorner1, nextradius, prevradius, threshold, curvedconvthreshold, curvednewtonfactor)
                            except:
                                if mergedorientations[path[i+1]]:
                                    print("Failed: getCurvedCornerFacet({}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(curvedcornerpolys, curvedcornerareas, curvedcorner1, curvedcornerpoint, curvedcorner2, prevradius, nextradius, threshold, curvedconvthreshold, curvednewtonfactor))
                                    curvedcorner = [curvedcorner1, curvedcornerpoint, curvedcorner2]
                                else:
                                    print("Failed: getCurvedCornerFacet({}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(curvedcornerpolys, curvedcornerareas, curvedcorner2, curvedcornerpoint, curvedcorner1, nextradius, prevradius, threshold, curvedconvthreshold, curvednewtonfactor))
                                    curvedcorner = [curvedcorner2, curvedcornerpoint, curvedcorner1]
                            for cornerindex in range(prevarc+1, nextarc):
                                if mergedorientations[path[i+1]]:
                                    facetfitted[cornerindex % len(facetfitted)] = ['curvedcorner', prevcenter, nextcenter, prevradius, nextradius, curvedcorner.copy()]
                                else:
                                    facetfitted[cornerindex % len(facetfitted)] = ['curvedcorner', nextcenter, prevcenter, nextradius, prevradius, curvedcorner.copy()]

                #Make facets gapless
                if makeGapless:
                    for i in range(len(facetfitted)):
                        if facetfitted[i] is not None and facetfitted[i][0] is not 'corner' and facetfitted[i][0] is not 'curvedcorner' and facetfitted[(i+1) % len(facetfitted)] is not None and facetfitted[(i+1) % len(facetfitted)][0] is not 'corner' and facetfitted[(i+1) % len(facetfitted)][0] is not 'curvedcorner':
                            #average rightmost intersect of current facet and leftmost intersect of previous facet
                            newIntersect = [(facetfitted[i][-1][-1][0]+facetfitted[(i+1) % len(facetfitted)][-1][0][0])/2, (facetfitted[i][-1][-1][1]+facetfitted[(i+1) % len(facetfitted)][-1][0][1])/2]
                            facetfitted[i][-1][-1] = newIntersect
                            facetfitted[(i+1)% len(facetfitted)][-1][0] = newIntersect
                    for i in range(len(facetfitted)):
                        #update circular facet center if needed
                        if facetfitted[i] is not None and facetfitted[i][0] == 'arc':
                            intersectleft = facetfitted[i][3][0]
                            intersectright = facetfitted[i][3][-1]
                            assert getDistance(intersectleft, intersectright) < 2*facetfitted[i][2]
                            #compute new center
                            #TODO: instead of preserving curvature, try preserving volume fraction
                            normal = [math.sqrt(facetfitted[i][2]**2 - getDistance(intersectleft, intersectright)**2/4)/getDistance(intersectleft, intersectright)*(intersectleft[1]-intersectright[1]), math.sqrt(facetfitted[i][2]**2 - getDistance(intersectleft, intersectright)**2/4)/getDistance(intersectleft, intersectright)*(intersectright[0]-intersectleft[0])]
                            facetcenter1 = [(intersectleft[0]+intersectright[0])/2 - normal[0], (intersectleft[1]+intersectright[1])/2 - normal[1]]
                            facetcenter2 = [(intersectleft[0]+intersectright[0])/2 + normal[0], (intersectleft[1]+intersectright[1])/2 + normal[1]]
                            if getDistance(facetfitted[i][1], facetcenter1) < getDistance(facetfitted[i][1], facetcenter2):
                                facetfitted[i][1] = facetcenter1
                            else:
                                facetfitted[i][1] = facetcenter2
                
                #Store facet info for each grid square (from merged squares)
                for i in range(len(facetfitted)):
                    for facetsquares in mergedpolyinfos[path[i+1]][0]:
                        predfacets[facetsquares[0]][facetsquares[1]] = facetfitted[i]

                #Get unique facets
                for facet in facetfitted:
                    if facet not in facetsunique:
                        facetsunique.append(facet)

    return predfacets, predorientations, facetsunique

def advectFacets(opolys, areas, predfacets, predorientations, velocity, threshold):
    #areas after new timestep

    totalunshiftedareas = []

    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if predfacets[x][y] is None and areas[x][y] > 0:
                totalunshiftedareas.append(areas[x][y]*getArea(opolys[x][y]))
            elif predfacets[x][y] is not None:
                predfacettype = predfacets[x][y][0]
                if predfacettype == 'linear':
                    totalunshiftedareas.append(getPolyLineArea(opolys[x][y], predfacets[x][y][1][0], predfacets[x][y][1][1]))
                elif predfacettype == 'corner':
                    cornerfacet = predfacets[x][y][1]
                    polyintersectioncornerarea = getPolyCornerArea(opolys[x][y], cornerfacet[0], cornerfacet[1], cornerfacet[2])*getArea(opolys[x][y])
                    if predorientations[x][y]:
                        totalunshiftedareas.append(abs(polyintersectioncornerarea))
                    else:
                        totalunshiftedareas.append(abs(getArea(opolys[x][y]))-abs(polyintersectioncornerarea))
                elif predfacettype == 'arc':
                    polyintersectioncirclearea, _ = getCircleIntersectArea([predfacets[x][y][1][0], predfacets[x][y][1][1]], predfacets[x][y][2], opolys[x][y])
                    if predorientations[x][y]:
                        totalunshiftedareas.append(abs(polyintersectioncirclearea))
                    else:
                        totalunshiftedareas.append(abs(getArea(opolys[x][y]))-abs(polyintersectioncirclearea))
                else: #curved corner:
                    curvedcorner = predfacets[x][y][-1]
                    polyintersectioncurvedcornerarea = getArea(opolys[x][y])*getPolyCurvedCornerArea(opolys[x][y], curvedcorner[0], curvedcorner[1], curvedcorner[2], predfacets[x][y][3], predfacets[x][y][4])
                    if predorientations[x][y]:
                        totalunshiftedareas.append(abs(polyintersectioncurvedcornerarea))
                    else:
                        totalunshiftedareas.append(abs(getArea(opolys[x][y]))-abs(polyintersectioncurvedcornerarea))
                #totalunshiftedareas[len(totalunshiftedareas)-1] /= getArea(opolys[x][y])

    print("Preshifted area: {}".format(sum(totalunshiftedareas)))

    totalgivenarea = 0
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            totalgivenarea += areas[x][y]*getArea(opolys[x][y])
    print("Given area: {}".format(totalgivenarea))

    totalshiftedareas = []

    nareas = [[[] for _ in range(len(opolys))] for _ in range(len(opolys))]
    
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if predfacets[x][y] is not None:
                predfacettype = predfacets[x][y][0]
                shiftpoly = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], opolys[x][y]))
                shiftbounds = [min(list(map(lambda x : x[0], shiftpoly))), min(list(map(lambda x : x[1], shiftpoly))), max(list(map(lambda x : x[0], shiftpoly))), max(list(map(lambda x : x[1], shiftpoly)))]
                #bounds: minx, miny, maxx, maxy
                #For each neighbor of the cell
                shiftareas = [[0]*3 for _ in range(3)]
                for testx in range(-1, 2):
                    for testy in range(-1, 2):
                        if x-testx >= 0 and x-testx < len(opolys) and y-testy >= 0 and y-testy < len(opolys):
                            testpoly = opolys[x-testx][y-testy]
                            testbounds = [min(list(map(lambda x : x[0], testpoly))), min(list(map(lambda x : x[1], testpoly))), max(list(map(lambda x : x[0], testpoly))), max(list(map(lambda x : x[1], testpoly)))]
                            if not(testbounds[2] <= shiftbounds[0] or shiftbounds[2] <= testbounds[0] or testbounds[3] <= shiftbounds[1] or shiftbounds[3] <= testbounds[1]):
                                #bounding boxes intersect, could be nonzero intersection
                                try:
                                    polyintersections = getPolyIntersectArea(shiftpoly, testpoly)
                                    if len(polyintersections) == 0:
                                        #No intersection
                                        continue
                                    #For each overlap region
                                    for polyintersection in polyintersections:
                                        if len(polyintersection) > 2:
                                            polyintersection = polyintersection[:len(polyintersection)-1]
                                            if getArea(polyintersection) < 0:
                                                polyintersection.reverse()
                                        else:
                                            continue
                                        
                                        if predfacettype == 'linear':
                                            shiftlinearfacet = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], predfacets[x][y][1]))
                                            polyintersectionlineararea = getPolyLineArea(polyintersection, shiftlinearfacet[0], shiftlinearfacet[len(shiftlinearfacet)-1])
                                            if predorientations[x][y]:
                                                shiftareas[testx+1][testy+1] += abs(polyintersectionlineararea)
                                            else:
                                                shiftareas[testx+1][testy+1] += abs(getArea(polyintersection))-abs(polyintersectionlineararea)
                                        elif predfacettype == 'corner':
                                            shiftcornerfacet = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], predfacets[x][y][-1]))
                                            polyintersectioncornerarea = getPolyCornerArea(polyintersection, shiftcornerfacet[0], shiftcornerfacet[1], shiftcornerfacet[2])*getArea(polyintersection)
                                            if predorientations[x][y]:
                                                shiftareas[testx+1][testy+1] += abs(polyintersectioncornerarea)
                                            else:
                                                shiftareas[testx+1][testy+1] += abs(getArea(polyintersection))-abs(polyintersectioncornerarea)
                                        elif predfacettype == 'arc':
                                            polyintersectioncirclearea, _ = getCircleIntersectArea([predfacets[x][y][1][0]+velocity[0], predfacets[x][y][1][1]+velocity[1]], predfacets[x][y][2], polyintersection)
                                            if predorientations[x][y]:
                                                shiftareas[testx+1][testy+1] += abs(polyintersectioncirclearea)
                                            else:
                                                shiftareas[testx+1][testy+1] += abs(getArea(polyintersection))-abs(polyintersectioncirclearea)
                                        else: #curved corner:
                                            curvedcorner = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], predfacets[x][y][-1]))
                                            polyintersectioncurvedcornerarea = getArea(polyintersection)*getPolyCurvedCornerArea(polyintersection, curvedcorner[0], curvedcorner[1], curvedcorner[2], predfacets[x][y][3], predfacets[x][y][4])
                                            if predorientations[x][y]:
                                                shiftareas[testx+1][testy+1] += abs(polyintersectioncurvedcornerarea)
                                            else:
                                                shiftareas[testx+1][testy+1] += abs(getArea(polyintersection))-abs(polyintersectioncurvedcornerarea)

                                except:
                                    print("Error in computing {} intersection of mesh cells: ({}, {}), ({}, {})".format(predfacettype, x, y, x-testx, y-testy))
                                    print("{}, {}".format(predfacets[x][y], polyintersection))

                for testx in range(len(shiftareas)):
                    for testy in range(len(shiftareas[0])):
                        nareas[x-testx+1][y-testy+1].append([shiftareas[testx][testy], [x, y]])

            elif areas[x][y] == 1:
                shiftpoly = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], opolys[x][y]))
                shiftbounds = [min(list(map(lambda x : x[0], shiftpoly))), min(list(map(lambda x : x[1], shiftpoly))), max(list(map(lambda x : x[0], shiftpoly))), max(list(map(lambda x : x[1], shiftpoly)))]
                shiftareas = [[0]*3 for _ in range(3)]
                #bounds: minx, miny, maxx, maxy
                for testx in range(-1, 2):
                    for testy in range(-1, 2):
                        if x-testx >= 0 and x-testx < len(opolys) and y-testy >= 0 and y-testy < len(opolys):
                            testpoly = opolys[x-testx][y-testy]
                            testbounds = [min(list(map(lambda x : x[0], testpoly))), min(list(map(lambda x : x[1], testpoly))), max(list(map(lambda x : x[0], testpoly))), max(list(map(lambda x : x[1], testpoly)))]
                            if not(testbounds[2] <= shiftbounds[0] or shiftbounds[2] <= testbounds[0] or testbounds[3] <= shiftbounds[1] or shiftbounds[3] <= testbounds[1]):
                                #bounding boxes intersect, could be nonzero intersection
                                try:
                                    polyintersections = getPolyIntersectArea(shiftpoly, testpoly)
                                    if len(polyintersections) == 0:
                                        #No intersection
                                        continue
                                    #assert len(polyintersections) == 1, print("getPolyIntersectArea({}, {})".format(shiftpoly, testpoly))
                                    for polyintersection in polyintersections:
                                        shiftareas[testx+1][testy+1] += abs(getArea(polyintersection))
                                except:
                                    print("Error in computing intersection of full mesh cells: ({}, {}), ({}, {})".format(x, y, x-testx, y-testy))
                                    
                for testx in range(len(shiftareas)):
                    for testy in range(len(shiftareas[0])):
                        if x-testx+1 >= 0 and x-testx+1 < len(opolys) and y-testy+1 >= 0 and y-testy+1 < len(opolys):
                            nareas[x-testx+1][y-testy+1].append([shiftareas[testx][testy], [x, y]])

            #Unpredicted partial fraction
            elif areas[x][y] > 0:
                print("Unpredicted partial fraction: {}, {}".format(x, y))
                
    """
    #plot errors in advecting facets
    plt.hist(list(map(lambda x : math.log10(x), totalshiftedareas)))
    plt.xlabel("|advected-original|")
    plt.ylabel("Frequency")
    plt.title("Error in advected cell area vs. original cell area")
    plt.savefig("advection_plots/advectarea_errors.png")
    plt.clf()
    """

    #sum immediate advected areas, could have issues
    newareas = [[0] * (len(opolys)) for _ in range(len(opolys))]
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            newareas[x][y] = sum(list(map(lambda z : z[0], nareas[x][y])))

    print("Sum after simple advection: {}".format(sum(list(map(lambda z : sum(z), newareas)))))

    #loop again through all advected areas, if any errors remain, fix them globally
    sumofexcesses = 0
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if newareas[x][y] > getArea(opolys[x][y]):
                sumofexcesses += newareas[x][y]-getArea(opolys[x][y])
                newareas[x][y] = getArea(opolys[x][y])
            elif newareas[x][y] != 0 and abs(newareas[x][y])/getArea(opolys[x][y]) < threshold:
                sumofexcesses += newareas[x][y]
                newareas[x][y] = 0
            elif newareas[x][y] != getArea(opolys[x][y]) and abs(1 - (newareas[x][y]/getArea(opolys[x][y]))) < threshold:
                sumofexcesses -= getArea(opolys[x][y])-newareas[x][y]
                newareas[x][y] = getArea(opolys[x][y])

    print("Sum of excesses: {}".format(sumofexcesses))
    
    print("Sum after fixing excesses: {}".format(sum(list(map(lambda x : sum(x), newareas)))))

    #Print errors in newareas
    """
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if not ((areas[x][y] == 0 and newareas[x][y]/getArea(opolys[x][y]) == 0) or (areas[x][y] == 1 and newareas[x][y]/getArea(opolys[x][y]) == 1)) and abs(areas[x][y] - newareas[x][y]/getArea(opolys[x][y])) > 1e-3:
                print("{}, {}, {}, {}".format([x, y], newareas[x][y]/getArea(opolys[x][y]), areas[x][y], abs(newareas[x][y]/getArea(opolys[x][y]) - areas[x][y])))
    """

    #Compute fraction to account for fixing excesses
    sumpartialgaps = 0
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if 0 < newareas[x][y] and newareas[x][y] < getArea(opolys[x][y]):
                sumpartialgaps += getArea(opolys[x][y])-newareas[x][y]
    fixexcessfraction = (totalgivenarea - sum(list(map(lambda x : sum(x), newareas))))/sumpartialgaps

    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if 0 < newareas[x][y] and newareas[x][y] < getArea(opolys[x][y]):
                newareas[x][y] = newareas[x][y]+(getArea(opolys[x][y])-newareas[x][y])*fixexcessfraction
            newareas[x][y] /= getArea(opolys[x][y])

    return newareas

    #--------------------------------

    #fix errors with area intersection computation
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            newarea = newareas[x][y]
            if newarea > getArea(opolys[x][y]):
                #print("Over 1: {}".format(newarea))
                #print(getArea(opolys[x][y]))
                #too much area in [x, y], move it somewhere else
                excessarea = newarea - getArea(opolys[x][y])
                advectionneighbors = []
                #advection neighbors with less than total fraction
                for potentialneighbor in nareas[x][y]:
                    if newareas[potentialneighbor[1][0]][potentialneighbor[1][1]] < getArea(opolys[potentialneighbor[1][0]][potentialneighbor[1][1]]) and newareas[potentialneighbor[1][0]][potentialneighbor[1][1]] > getArea(opolys[potentialneighbor[1][0]][potentialneighbor[1][1]])*threshold and potentialneighbor[0] > 0:
                        advectionneighbors.append(potentialneighbor)
                while excessarea > 0 and len(advectionneighbors) > 0:
                    #repeatedly remove excess area from [x, y] in proportions of advected areas until no excess area or no neighbors
                    minadvectionfraction = float('inf')
                    minadvectionneighbor = None
                    advectionneighborareasum = sum(list(map(lambda z : z[0], advectionneighbors)))
                    for advectionneighborinfo in advectionneighbors:
                        advectionneighbor = advectionneighborinfo[1]
                        advectionneighborareagap = getArea(opolys[advectionneighbor[0]][advectionneighbor[1]]) - newareas[advectionneighbor[0]][advectionneighbor[1]]
                        currentadvectionfraction = advectionneighborareagap/(excessarea*(advectionneighborinfo[0])/advectionneighborareasum)
                        if currentadvectionfraction < minadvectionfraction:
                            minadvectionfraction = currentadvectionfraction
                            minadvectionneighbor = advectionneighborinfo
                    if minadvectionfraction > 1:
                        minadvectionfraction = 1
                        minadvectionneighbor = None
                    for advectionneighborinfo in advectionneighbors:
                        advectionneighbor = advectionneighborinfo[1]
                        newareas[advectionneighbor[0]][advectionneighbor[1]] += minadvectionfraction*excessarea*(advectionneighborinfo[0])/advectionneighborareasum
                    newareas[x][y] -= excessarea * minadvectionfraction
                    excessarea *= 1-minadvectionfraction
                    if minadvectionneighbor is not None:
                        advectionneighbors.remove(minadvectionneighbor)
                #after the loop, either no more excess area or all neighbors are full
                """
                if excessarea > 0:
                    sumofexcesses += excessarea
                newareas[x][y] = 1
                """
                #print(newareas[x][y])
            elif newarea != 0 and abs(newarea)/getArea(opolys[x][y]) < threshold:
                #print("Almost 0: {}".format(newarea))
                #print(0)
                #treating [x, y] as 0 area, move area somewhere else
                excessarea = newarea
                advectionneighbors = []
                #advection neighbors with less than total fraction
                for potentialneighbor in nareas[x][y]:
                    if newareas[potentialneighbor[1][0]][potentialneighbor[1][1]] < getArea(opolys[potentialneighbor[1][0]][potentialneighbor[1][1]]) and newareas[potentialneighbor[1][0]][potentialneighbor[1][1]] > getArea(opolys[potentialneighbor[1][0]][potentialneighbor[1][1]])*threshold and potentialneighbor[0] > 0:
                        advectionneighbors.append(potentialneighbor)
                while excessarea > 0 and len(advectionneighbors) > 0:
                    #repeatedly remove excess area from [x, y] in proportions of advected areas until no excess area or no neighbors
                    minadvectionfraction = float('inf')
                    minadvectionneighbor = None
                    advectionneighborareasum = sum(list(map(lambda z : z[0], advectionneighbors)))
                    for advectionneighborinfo in advectionneighbors:
                        advectionneighbor = advectionneighborinfo[1]
                        advectionneighborareagap = getArea(opolys[advectionneighbor[0]][advectionneighbor[1]]) - newareas[advectionneighbor[0]][advectionneighbor[1]]
                        currentadvectionfraction = advectionneighborareagap/(excessarea*(advectionneighborinfo[0])/advectionneighborareasum)
                        if currentadvectionfraction < minadvectionfraction:
                            minadvectionfraction = currentadvectionfraction
                            minadvectionneighbor = advectionneighborinfo
                    if minadvectionfraction > 1:
                        minadvectionfraction = 1
                        minadvectionneighbor = None
                    for advectionneighborinfo in advectionneighbors:
                        advectionneighbor = advectionneighborinfo[1]
                        newareas[advectionneighbor[0]][advectionneighbor[1]] += minadvectionfraction*excessarea*(advectionneighborinfo[0])/advectionneighborareasum
                    newareas[x][y] -= excessarea * minadvectionfraction
                    excessarea *= 1-minadvectionfraction
                    if minadvectionneighbor is not None:
                        advectionneighbors.remove(minadvectionneighbor)
                #after the loop, either no more excess area or all neighbors are full
                """
                if excessarea > 0:
                    sumofexcesses += excessarea
                newareas[x][y] = 0
                """
                #print(newareas[x][y])
            elif newarea != getArea(opolys[x][y]) and abs((newarea/getArea(opolys[x][y])) - 1) < threshold:
                #print("Almost 1: {}".format(newarea))
                #print(getArea(opolys[x][y]))
                #treating [x, y] as 1 area, needs to take area from somewhere else
                neededarea = 1-newarea
                advectionneighbors = []
                #advection neighbors with more than total fraction
                for potentialneighbor in nareas[x][y]:
                    if newareas[potentialneighbor[1][0]][potentialneighbor[1][1]] > getArea(opolys[potentialneighbor[1][0]][potentialneighbor[1][1]]) and newareas[potentialneighbor[1][0]][potentialneighbor[1][1]] < getArea(opolys[potentialneighbor[1][0]][potentialneighbor[1][1]])*(1-threshold) and potentialneighbor[0] > 0:
                        advectionneighbors.append(potentialneighbor)
                while neededarea > 0 and len(advectionneighbors) > 0:
                    #repeatedly remove excess area from [x, y] in proportions of advected areas until no excess area or no neighbors
                    minadvectionfraction = float('inf')
                    minadvectionneighbor = None
                    advectionneighborareasum = sum(list(map(lambda z : z[0], advectionneighbors)))
                    for advectionneighborinfo in advectionneighbors:
                        advectionneighbor = advectionneighborinfo[1]
                        advectionneighborareagap = newareas[advectionneighbor[0]][advectionneighbor[1]] - getArea(opolys[advectionneighbor[0]][advectionneighbor[1]])
                        currentadvectionfraction = advectionneighborareagap/(neededarea*(advectionneighborinfo[0])/advectionneighborareasum)
                        if currentadvectionfraction < minadvectionfraction:
                            minadvectionfraction = currentadvectionfraction
                            minadvectionneighbor = advectionneighborinfo
                    if minadvectionfraction > 1:
                        minadvectionfraction = 1
                        minadvectionneighbor = None
                    for advectionneighborinfo in advectionneighbors:
                        advectionneighbor = advectionneighborinfo[1]
                        newareas[advectionneighbor[0]][advectionneighbor[1]] -= minadvectionfraction*neededarea*(advectionneighborinfo[0])/advectionneighborareasum
                    newareas[x][y] += neededarea * minadvectionfraction
                    neededarea *= 1-minadvectionfraction
                    if minadvectionneighbor is not None:
                        advectionneighbors.remove(minadvectionneighbor)
                #after the loop, either no more excess area or all neighbors are full
                """
                if excessarea > 0:
                    sumofexcesses += excessarea
                newareas[x][y] = 0
                """
                #print(newareas[x][y])
            #print(sum(list(map(lambda x : sum(x), newareas))))

    print("Sum before fixing excesses: {}".format(sum(list(map(lambda x : sum(x), newareas)))))

    #loop again through all advected areas, if any errors remain, fix them globally
    sumofexcesses = 0
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if newareas[x][y] > getArea(opolys[x][y]):
                sumofexcesses += newareas[x][y]-getArea(opolys[x][y])
                newareas[x][y] = getArea(opolys[x][y])
            elif newareas[x][y] != 0 and abs(newareas[x][y])/getArea(opolys[x][y]) < threshold:
                sumofexcesses += newareas[x][y]
                newareas[x][y] = 0
            elif newareas[x][y] != getArea(opolys[x][y]) and abs(1 - (newareas[x][y]/getArea(opolys[x][y]))) < threshold:
                sumofexcesses -= getArea(opolys[x][y])-newareas[x][y]
                newareas[x][y] = getArea(opolys[x][y])

    print("Sum of excesses: {}".format(sumofexcesses))
    
    print("Sum after fixing excesses: {}".format(sum(list(map(lambda x : sum(x), newareas)))))

    #Compute fraction to account for fixing excesses
    sumpartialgaps = 0
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if 0 < newareas[x][y] and newareas[x][y] < getArea(opolys[x][y]):
                sumpartialgaps += getArea(opolys[x][y])-newareas[x][y]
    fixexcessfraction = (totalgivenarea - sum(list(map(lambda x : sum(x), newareas))))/sumpartialgaps

    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if 0 < newareas[x][y] and newareas[x][y] < getArea(opolys[x][y]):
                newareas[x][y] = newareas[x][y]+(getArea(opolys[x][y])-newareas[x][y])*fixexcessfraction
            newareas[x][y] /= getArea(opolys[x][y])

    #print("Sum after fixing areas after excesses: {}".format(sum(list(map(lambda x : sum(x), newareas)))))

    return newareas

def advection(gridSize, timesteps, plotVtk, plotMat, makeGapless):
    
    #starting mesh
    #opoints = makeCartesianGrid(gridSize)
    #opoints = makeQuadGrid(gridSize, wiggle)
    #opoints = makeConcaveGrid(gridSize, wiggle)
    opoints = makeFineCartesianGrid(gridSize, resolution)
    
    #mesh velocity per timestep
    velocity = [vsize, vsize]
    
    #variables per time step
    opolys = [[0] * (len(opoints)-1) for _ in range(len(opoints)-1)]

    #matplotlib
    opatches = []
    opatchareas = []
    opatchareapartials = []
    
    #compute initial area setup
    areas = [[0] * len(opolys) for _ in range(len(opolys))]

    for x in range(len(opoints)-1):
        for y in range(len(opoints)-1):
            opoly = [opoints[x][y], opoints[x+1][y], opoints[x+1][y+1], opoints[x][y+1]]
            
            #Initial area setting
            """
            #cross
            xpoints = [[4, 0], [10, 6], [16, 0], [20, 4], [14, 10], [20, 16], [16, 20], [10, 14], [4, 20], [0, 16], [6, 10], [0, 4]]
            xpoints = list(map(lambda x : [x[0]+3.005, x[1]+3.025], xpoints))
            xpolyintersects = getPolyIntersectArea(xpoints, opoly)
            for xpolyintersect in xpolyintersects:
                areas[x][y] += abs(getArea(xpolyintersect))
            #cross
            xpoints = [[6, 0], [14, 0], [14, 6], [20, 6], [20, 14], [14, 14], [14, 20], [6, 20], [6, 14], [0, 14], [0, 6], [6, 6]]
            xpoints = list(map(lambda x : [x[0]+3.005, x[1]+25.025], xpoints))
            xpolyintersects = getPolyIntersectArea(xpoints, opoly)
            for xpolyintersect in xpolyintersects:
                areas[x][y] += abs(getArea(xpolyintersect))
            
            
            truecircleintersect1, truecircleintersect2 = getCircleCircleIntersects([38.005, 33.005], [38.005, 13.005], 12, 12)
            if getDistance(getCentroid(opoly), truecircleintersect1) <= 1:
                areas[x][y] += getPolyCurvedCornerArea(opoly, [26.005, 33.005], truecircleintersect1, [26.005, 13.005], 12, 12)*getArea(opoly)
            elif getDistance(getCentroid(opoly), truecircleintersect1) <= 1:#(x >= 43 and x <= 45 and y >= 21 and y <= 23):
                areas[x][y] += getPolyCurvedCornerArea(opoly, [50.005, 13.005], truecircleintersect2, [50.005, 33.005], 12, 12)*getArea(opoly)
                print("{}, {}".format(opoly, getPolyCurvedCornerArea(opoly, [50.005, 13.005], truecircleintersect2, [50.005, 33.005], 12, 12)))
            elif getDistance(getCentroid(opoly), [38.005, 33.005]) < getDistance(getCentroid(opoly), [38.005, 13.005]):
                center = [38.005, 33.005]
                area, intersect = getCircleIntersectArea(center, 12, opoly)
                areas[x][y] += area
            else:
                center = [38.005, 13.005]
                area, intersect = getCircleIntersectArea(center, 12, opoly)
                areas[x][y] += area
            
            #ring
            radius = 12
            radiussmall = 7
            center = [38.005, 13.005]
            #area, intersect = getCircleIntersectArea(center, radius, opoly)
            #areas[x][y] += area
            area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
            areas[x][y] -= area
            #ring
            radius = 12
            radiussmall = 7
            center = [38.005, 33.005]
            #area, intersect = getCircleIntersectArea(center, radius, opoly)
            #areas[x][y] += area
            area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
            areas[x][y] -= area
            """

            #rings and dot
            """
            #ring
            radius = 12
            radiussmall = 7
            center = [38.005, 13.005]
            area, intersect = getCircleIntersectArea(center, radius, opoly)
            areas[x][y] += area
            area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
            areas[x][y] -= area
            #ring
            radius = 12
            radiussmall = 7
            center = [38.005, 43.005]
            area, intersect = getCircleIntersectArea(center, radius, opoly)
            areas[x][y] += area
            area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
            areas[x][y] -= area

            #dot
            radius = 3
            center = [38.005, 13.005]
            area, intersect = getCircleIntersectArea(center, radius, opoly)
            areas[x][y] += area
            """

            #ellipse
            #"""
            ellipsecenter = [16.000005, 16.000005]
            ellipseradius = 7
            majortimes = 1.5
            shiftedpoly = list(map(lambda x : [x[0]/majortimes, x[1]], opoly))
            area, intersect = getCircleIntersectArea([ellipsecenter[0]/majortimes, ellipsecenter[1]], ellipseradius, shiftedpoly)
            areas[x][y] += majortimes*area
            #"""

            areas[x][y] /= getArea(opoly)
            if areas[x][y] > 1:
                areas[x][y] = 1
            if abs(1-areas[x][y]) < threshold:
                areas[x][y] = 1
            elif abs(areas[x][y]) < threshold:
                areas[x][y] = 0
            elif areas[x][y] < 0:
                areas[x][y] = 0

            opolys[x][y] = opoly

            if plotMat:
                #matplotlib: match partial areas and polygons for plotting
                opatch = Polygon(np.array(opoly), True)
                opatches.append(opatch)
                opatchareas.append(areas[x][y])
                opatchareapartials.append(math.ceil(areas[x][y] - math.floor(areas[x][y])))

    #plot in matplotlib
    if plotMat:
        try:
            os.mkdir('advection_plots')
        except:
            print("Saving plots in ./advection_plots/.")
        op = PatchCollection(opatches, cmap='jet')
        opatchareas = np.array(opatchareas)
        op.set_array(opatchareas)
        fig, ax = plt.subplots()
        ax.set_xlim(-1, len(opolys)/resolution+1)
        ax.set_ylim(-1, len(opolys)/resolution+1)
        ax.add_collection(op)
        plt.savefig("advection_plots/original.png", dpi=199)
        plt.clf()
        
        op = PatchCollection(opatches, cmap='jet')
        opatchareas = np.array(opatchareapartials)
        op.set_array(opatchareas)
        fig, ax = plt.subplots()
        ax.set_xlim(-1, len(opolys)/resolution+1)
        ax.set_ylim(-1, len(opolys)/resolution+1)
        ax.add_collection(op)
        plt.savefig("advection_plots/original_partials.png", dpi=199)
        plt.clf()
    
    """
    #plot in vtk
    if plotVtk:
        try:
            os.mkdir('advection_vtk')
        except:
            print("Saving vtk files in ./advection_vtk/.")
        plotQuadGrid(opoints, 'advection_vtk/quads')
    """
    
    #main advection loop
    for timestep in range(timesteps):
        print("Timestep: {}".format(timestep))
        
        predmergedpolys, mergedpolys, mergedcoords, mergedareafractions = merge(opolys, areas)
        
        predfacets, predorientations, facetsunique = makeFacets(predmergedpolys, mergedpolys, mergedcoords, mergedareafractions, opolys, areas, threshold, makeGapless)

        #curvature
        """
        curvatures1 = []
        curvatures2 = []
        for predfacet in facetsunique:
            try:
                facetintersect1 = predfacet[3][0]
                facetintersect2 = predfacet[3][-1]
                #facetintersect = [(facetintersect1[0]+facetintersect2[0])/2, (facetintersect1[1]+facetintersect2[1])/2]
                ellipsex = (facetintersect1[0] - ellipsecenter[0])/getDistance(facetintersect1, [ellipsecenter[0], ellipsecenter[1]])
                ellipsey = (facetintersect1[1] - ellipsecenter[1])/getDistance(facetintersect1, [ellipsecenter[0], ellipsecenter[1]])
                ellipsecurvature1 = ellipseradius**2 * majortimes/math.sqrt((ellipseradius*majortimes*ellipsey)**2 + (ellipseradius*ellipsex)**2)**3
                facetcurvature = 1/predfacet[2]
                curvatures1.append(max(1e-15, abs(ellipsecurvature1-facetcurvature)))
                ellipsex = (facetintersect2[0] - ellipsecenter[0])/getDistance(facetintersect2, [ellipsecenter[0], ellipsecenter[1]])
                ellipsey = (facetintersect2[1] - ellipsecenter[1])/getDistance(facetintersect2, [ellipsecenter[0], ellipsecenter[1]])
                ellipsecurvature2 = ellipseradius**2 * majortimes/math.sqrt((ellipseradius*majortimes*ellipsey)**2 + (ellipseradius*ellipsex)**2)**3
                curvatures2.append(max(1e-15, abs(ellipsecurvature2-facetcurvature)))
                print("{}, {}, {}".format(ellipsecurvature1, facetcurvature, ellipsecurvature2))
            except:
                pass
        print(curvatures1)
        plt.hist(list(map(lambda x : math.log10(x), curvatures1)))
        plt.xlabel("Log error in curvature")
        plt.ylabel("Frequency")
        plt.title("Error in curvature")
        plt.savefig("advection_plots/ellipse_curvature1_errors.png")
        plt.clf()
        print(curvatures2)
        plt.hist(list(map(lambda x : math.log10(x), curvatures2)))
        plt.xlabel("Log error in curvature")
        plt.ylabel("Frequency")
        plt.title("Error in curvature")
        plt.savefig("advection_plots/ellipse_curvature2_errors.png")
        plt.clf()
        """

        #Plot facets in vtk
        if plotVtk and (timestep % 1 == 0 or timestep == timesteps-1):
            plotFacets(facetsunique, 'advection_vtk/timestep_{}'.format(timestep))
        
        print("Computing new areas")
        nareas = advectFacets(opolys, areas, predfacets, predorientations, velocity, threshold)
        areas = nareas
        
        """
        prederrors = []
        for x in range(len(areas)):
            for y in range(len(areas)):
                center7 = [38.005+0.09*(timestep+1), 13.005+0.09*(timestep+1)]
                center8 = [38.005+0.09*(timestep+1), 33.005+0.09*(timestep+1)]
                if getDistance([x+0.5, y+0.5], center7) < getDistance([x+0.5, y+0.5], center8):
                    testarea9, _ = getCircleIntersectArea(center7, 12, [[x, y], [x+1, y], [x+1, y+1], [x, y+1]])
                else:
                    testarea9, _ = getCircleIntersectArea(center8, 12, [[x, y], [x+1, y], [x+1, y+1], [x, y+1]])
                if abs(testarea9-areas[x][y]) > 1e-9:
                    print("{}, {}".format([x, y], abs(testarea9-areas[x][y])))
        
        plt.hist(list(map(lambda x : math.log10(x), prederrors)))
        plt.xlabel("Predicted area log error")
        plt.ylabel("Frequency")
        plt.title("Errors in predicted area, timestep={}".format(timestep))
        plt.savefig("advection_plots/area_errors_timestep_{}.png".format(timestep))
        plt.clf()
        """

        if plotMat:
            print("Plotting")
            #plot in matplotlib
            newp = PatchCollection(opatches, cmap='jet')
            npatchareas = []
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    npatchareas.append(areas[x][y])
            npatchareas = np.array(npatchareas)
            newp.set_array(npatchareas)
            fig, ax = plt.subplots()
            ax.set_xlim(-1, len(opolys)/resolution+1)
            ax.set_ylim(-1, len(opolys)/resolution+1)
            ax.add_collection(newp)
            plt.savefig("advection_plots/pred_timestep_{}.png".format(timestep), dpi=199)
            plt.clf()

            newp = PatchCollection(opatches, cmap='jet')
            npatchareas = []
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    npatchareas.append(math.ceil(areas[x][y] - math.floor(areas[x][y])))
            npatchareas = np.array(npatchareas)
            newp.set_array(npatchareas)
            fig, ax = plt.subplots()
            ax.set_xlim(-1, len(opolys)/resolution+1)
            ax.set_ylim(-1, len(opolys)/resolution+1)
            ax.add_collection(newp)
            plt.savefig("advection_plots/predpartials_timestep_{}.png".format(timestep), dpi=199)
            plt.clf()

            
            #True areas
            trueareas = [[0] * len(opolys) for _ in range(len(opolys))]
            truepatchareas = []
            prederrors = []
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    opoly = [opoints[x][y], opoints[x+1][y], opoints[x+1][y+1], opoints[x][y+1]]
                    """
                    #material 1 cross
                    xpoints = [[4, 0], [10, 6], [16, 0], [20, 4], [14, 10], [20, 16], [16, 20], [10, 14], [4, 20], [0, 16], [6, 10], [0, 4]]
                    xpoints = list(map(lambda x : [x[0]+3.005+velocity[0]*(timestep+1), x[1]+3.025+velocity[1]*(timestep+1)], xpoints))
                    xpolyintersects = getPolyIntersectArea(xpoints, opoly)
                    for xpolyintersect in xpolyintersects:
                        trueareas[x][y] += abs(getArea(xpolyintersect))
                    #material 2 cross
                    xpoints = [[6, 0], [14, 0], [14, 6], [20, 6], [20, 14], [14, 14], [14, 20], [6, 20], [6, 14], [0, 14], [0, 6], [6, 6]]
                    xpoints = list(map(lambda x : [x[0]+3.005+velocity[0]*(timestep+1), x[1]+25.025+velocity[1]*(timestep+1)], xpoints))
                    xpolyintersects = getPolyIntersectArea(xpoints, opoly)
                    for xpolyintersect in xpolyintersects:
                        trueareas[x][y] += abs(getArea(xpolyintersect))
                    
                    truecircleintersect1, truecircleintersect2 = getCircleCircleIntersects([38.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)], [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)], 12, 12)
                    if getDistance([x+0.5, y+0.5], truecircleintersect1) <= 2:
                        trueareas[x][y] += getPolyCurvedCornerArea(opoly, [26.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)], truecircleintersect1, [26.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)], 12, 12)
                    elif getDistance([x+0.5, y+0.5], truecircleintersect2) <= 2:
                        trueareas[x][y] += getPolyCurvedCornerArea(opoly, [50.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)], truecircleintersect2, [50.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)], 12, 12)
                    elif getDistance([x+0.5, y+0.5], [38.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)]) < getDistance([x+0.5, y+0.5], [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)]):
                        center = [38.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)]
                        area, intersect = getCircleIntersectArea(center, 12, opoly)
                        trueareas[x][y] += area
                    else:
                        center = [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)]
                        area, intersect = getCircleIntersectArea(center, 12, opoly)
                        trueareas[x][y] += area
                    
                    #material 1 ring
                    radius = 12
                    radiussmall = 7
                    center = [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)]
                    #area, intersect = getCircleIntersectArea(center, radius, opoly)
                    #trueareas[x][y] += area
                    area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
                    trueareas[x][y] -= area
                    #material 2 ring
                    radius = 12
                    radiussmall = 7
                    center = [38.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)]
                    #area, intersect = getCircleIntersectArea(center, radius, opoly)
                    #trueareas[x][y] += area
                    area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
                    trueareas[x][y] -= area
                    """

                    """
                    #material 1 ring
                    radius = 12
                    radiussmall = 7
                    center = [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)]
                    area, intersect = getCircleIntersectArea(center, radius, opoly)
                    trueareas[x][y] += area
                    area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
                    trueareas[x][y] -= area
                    #material 2 ring
                    radius = 12
                    radiussmall = 7
                    center = [38.005+velocity[0]*(timestep+1), 43.005+velocity[1]*(timestep+1)]
                    area, intersect = getCircleIntersectArea(center, radius, opoly)
                    trueareas[x][y] += area
                    area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
                    trueareas[x][y] -= area

                    #material 2 dot
                    radius = 3
                    center = [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)]
                    area, intersect = getCircleIntersectArea(center, radius, opoly)
                    trueareas[x][y] += area
                    """

                    ellipsecenter = [16.000005+velocity[0]*(timestep+1), 16.000005+velocity[1]*(timestep+1)]
                    shiftedpoly = list(map(lambda x : [x[0]/majortimes, x[1]], opoly))
                    area, intersect = getCircleIntersectArea([ellipsecenter[0]/majortimes, ellipsecenter[1]], ellipseradius, shiftedpoly)
                    trueareas[x][y] += majortimes*area

                    trueareas[x][y] /= getArea(opolys[x][y])
                    if abs(1-trueareas[x][y]) < threshold:
                        trueareas[x][y] = 1
                    elif abs(trueareas[x][y]) < threshold:
                        trueareas[x][y] = 0
                    truepatchareas.append(trueareas[x][y])
                    
                    if abs(trueareas[x][y] - areas[x][y]) > 0:
                        prederrors.append(abs(trueareas[x][y] - areas[x][y]))
                        if abs(trueareas[x][y] - areas[x][y]) > 1e-10:
                            print(abs(trueareas[x][y] - areas[x][y]))
            
            newp = PatchCollection(opatches, cmap='jet')
            npatchareas = np.array(truepatchareas)
            newp.set_array(npatchareas)
            fig, ax = plt.subplots()
            ax.set_xlim(-1, len(opolys)/resolution+1)
            ax.set_ylim(-1, len(opolys)/resolution+1)
            ax.add_collection(newp)
            plt.savefig("advection_plots/true_timestep_{}.png".format(timestep), dpi=199)
            plt.clf()
            #"""
            plt.hist(list(map(lambda x : math.log10(x), prederrors)))
            plt.xlabel("Predicted area fraction log error")
            plt.ylabel("Frequency")
            plt.title("Errors in predicted area fraction, timestep={}".format(timestep))
            plt.savefig("advection_plots/area_errors_timestep_{}.png".format(timestep))
            plt.clf()
            #"""

            prederrors = list(map(lambda x : max(1e-15, x), prederrors))
            prederrors = list(map(lambda x : math.log10(x), prederrors))
            print(sum(prederrors)/len(prederrors))
            totaladvectedarea = 0
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    totaladvectedarea += areas[x][y]*getArea(opolys[x][y])

            totaltruearea = 0
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    totaltruearea += trueareas[x][y]*getArea(opolys[x][y])
                    
            print("Advected area: {}".format(totaladvectedarea))
            print("True area: {}".format(totaltruearea))

            plt.close("all")

np.random.seed(17)
wiggle = 0.5
threshold = 1e-10
vsize = 0.0905
linearerrorthreshold = 0#1e-8
cornererrorthreshold = 1e-7
similararcthreshold = 1e-3
curvedconvthreshold = 1e-13
curvednewtonfactor = 0.1
gridSize = 50
resolution = 1
timesteps = 1
advection(gridSize, timesteps, plotVtk=False, plotMat=True, makeGapless=True)
#advection(gridSize, timesteps, plotVtk=True, plotMat=True, makeGapless=True)