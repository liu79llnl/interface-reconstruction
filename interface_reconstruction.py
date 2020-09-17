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
                    mergesquares.append([curx, cury])
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

#Hyperparameters:
linearerrorthreshold = 1e-8 #if area fraction error in linear facet < linearerrorthreshold, use linear facet at this cell
cornererrorthreshold = 1e-7 #if area fraction error in corner facet < cornererrorthreshold, use (straight edged) corner facet at this cell
curvaturethreshold = 1e-3 #if adjacent curvatures > curvaturethreshold, try to fit a curved corner facet
curvedcornerthreshold = 1e-3 #if area fraction error in curved corner < curvedcornerthreshold, use curved corner at this cell
curvedconvthreshold = 1e-10 #if curved corner optimization is only improving area fractions by < curvedconvthreshold, declare it converged
curvednewtonfactor = 0.1 #how much to step in Newton direction when optimizing curved corners

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
                        raise Exception("Interface path is unclear: increase fineness of resolution!")
                        print(mergedpolyinfos[curmergedpolyindex])
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
                
                if not(curOrientation):
                    path.reverse()
                    
                facetfitted = [None] * (len(path)-2)
                
                #Make linear facets
                for pathelement in range(1, len(path)-1):
                    #previ, curi, nexti are indices of mergedpolyinfos
                    curi = path[pathelement]
                    previ = path[pathelement-1]
                    nexti = path[pathelement+1]
                    prevpolygon = mergedcoords[previ]
                    prevpolygonarea = mergedareafractions[previ]
                    nextpolygon = mergedcoords[nexti]
                    nextpolygonarea = mergedareafractions[nexti]
                    try:
                        facetline1, facetline2 = getLinearFacet(prevpolygon, nextpolygon, prevpolygonarea, nextpolygonarea, threshold)
                        if abs(mergedareafractions[curi] - getPolyLineArea(mergedcoords[curi], facetline1, facetline2)/getArea(mergedcoords[curi])) > linearerrorthreshold:
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
                        prevFacet = lefts[pathelement]
                        nextFacet = rights[pathelement]
                        corner, _, _ = lineIntersect(prevFacet[0], prevFacet[1], nextFacet[0], nextFacet[1])
                        corner = [prevFacet[0], corner, nextFacet[1]]
                        cornerareafraction = getPolyCornerArea(mergedcoords[path[pathelement]], corner[0], corner[1], corner[2])/getArea(mergedcoords[path[pathelement]])
                        if abs(cornerareafraction - mergedareafractions[path[pathelement]]) < cornererrorthreshold:
                            facetfitted[pathelement] = ['corner', corner]
                        else:
                            continue

                #Try to make arc facets for the remaining ones
                for pathelement in range(1, len(path)-1):
                    if facetfitted[pathelement-1] is None:
                        #previ, curi, nexti are indices of mergedpolyinfos
                        curi = path[pathelement]
                        previ = path[pathelement-1]
                        nexti = path[pathelement+1]
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
                        previ = path[i]
                        nexti = path[i+2]
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
                for i in range(len(facetfitted)):
                    curi = path[i+1]
                    previ = path[i]
                    nexti = path[i+2]
                    curpolygon = mergedcoords[curi]
                    if facetfitted[i] is not None and facetfitted[i][0] == 'arc' and not(facetfitted[(i-1) % len(facetfitted)] is not None and (facetfitted[(i-1) % len(facetfitted)][0] == 'curvedcorner' or (facetfitted[(i-1) % len(facetfitted)][0] == 'arc' and abs((1/facetfitted[(i-1) % len(facetfitted)][2]) - (1/facetfitted[i][2])) < curvaturethreshold)) and facetfitted[(i+1) % len(facetfitted)] is not None and (facetfitted[(i+1) % len(facetfitted)][0] == 'curvedcorner' or (facetfitted[(i+1) % len(facetfitted)][0] == 'arc' and abs((1/facetfitted[(i+1) % len(facetfitted)][2]) - (1/facetfitted[i][2])) < curvaturethreshold))):
                        #Local curvatures differ by more than curvaturethreshold: try a curved corner
                        #Find leftward arc neighbor
                        prevarc = (i-1) % len(facetfitted)
                        #Go until you find two arc neighboring arc facets that have similar curvature
                        while prevarc != i and (facetfitted[prevarc] is None or facetfitted[(prevarc-1) % len(facetfitted)] is None or (facetfitted[prevarc][0] != 'arc' and facetfitted[prevarc][0] != 'linear') or (facetfitted[prevarc][0] == 'arc' and facetfitted[(prevarc-1) % len(facetfitted)][0] == 'arc' and abs(facetfitted[prevarc][2]-facetfitted[(prevarc-1) % len(facetfitted)][2]) > curvaturethreshold)):
                            prevarc = (prevarc-1) % len(facetfitted)
                        if prevarc == i:
                            #No neighbors
                            continue
                        nextarc = prevarc+1
                        #Go until you find two arc neighboring arc facets that have similar curvature
                        while nextarc != prevarc+len(facetfitted) and (facetfitted[nextarc % len(facetfitted)] is None or facetfitted[(nextarc+1) % len(facetfitted)] is None or (facetfitted[nextarc % len(facetfitted)][0] != 'arc' and facetfitted[nextarc % len(facetfitted)][0] != 'linear') or (facetfitted[nextarc % len(facetfitted)][0] == 'arc' and facetfitted[(nextarc+1) % len(facetfitted)][0] == 'arc' and abs(facetfitted[nextarc % len(facetfitted)][2]-facetfitted[(nextarc+1) % len(facetfitted)][2]) > curvaturethreshold)):
                            nextarc += 1
                        if nextarc % len(facetfitted) == prevarc:
                            #Only one neighbor
                            continue
                        nextarc = nextarc % len(facetfitted)
                        #Widen cluster of curved corner polys by 1 on both sides
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
                            curvedcornerpoint, _, _ = lineIntersect(prevcornerfacetside[1][0], prevcornerfacetside[1][1], nextcornerfacetside[1][0], nextcornerfacetside[1][1])
                        
                        #Use all polys in between two nearest successful facets
                        curvedcornergridindices = list(map(lambda x : path[(x % len(facetfitted))+1], range(prevarc+1, nextarc)))
                        curvedcornerpolys = list(map(lambda x : mergedcoords[x], curvedcornergridindices))
                        curvedcornerareas = list(map(lambda x : mergedareafractions[x], curvedcornergridindices))

                        curvedcornerguessareas = list(map(lambda x : getPolyCurvedCornerArea(mergedcoords[x], curvedcorner1, curvedcornerpoint, curvedcorner2, prevradius, nextradius)/getArea(mergedcoords[x]), curvedcornergridindices))
                        if max(list(map(lambda x : abs(curvedcornerareas[x]-curvedcornerguessareas[x]), range(len(curvedcornergridindices))))) < curvedcornerthreshold:
                            #Initial guess is close enough: optimize with Newton's method
                            try:
                                #curvedcorner = [curvedcorner1, curvedcornerpoint, curvedcorner2]
                                curvedcorner = getCurvedCornerFacet(curvedcornerpolys, curvedcornerareas, curvedcorner1, curvedcornerpoint, curvedcorner2, prevradius, nextradius, threshold, curvedconvthreshold, curvednewtonfactor)
                            except:
                                print("Failed: getCurvedCornerFacet({}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(curvedcornerpolys, curvedcornerareas, curvedcorner1, curvedcornerpoint, curvedcorner2, prevradius, nextradius, threshold, curvedconvthreshold, curvednewtonfactor))
                                curvedcorner = [curvedcorner1, curvedcornerpoint, curvedcorner2]
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
                            if facetfitted[i][2] > 0:
                                facetfitted[i][1] = getCenter(intersectleft, intersectright, facetfitted[i][2])
                            else:
                                facetfitted[i][1] = getCenter(intersectright, intersectleft, -facetfitted[i][2])
                            #compute new center
                            #TODO: instead of preserving curvature, try preserving volume fraction

                #Store facet info for each grid square (from merged squares)
                for i in range(len(facetfitted)):
                    for facetsquares in mergedpolyinfos[path[i+1]][0]:
                        predfacets[facetsquares[0]][facetsquares[1]] = facetfitted[i]

                #Get unique facets
                for facet in facetfitted:
                    if facet not in facetsunique:
                        facetsunique.append(facet)

    #print(facetsunique)

    return predfacets, facetsunique

def advectFacets(opolys, areas, predfacets, velocity, checksize, threshold):
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
                    polyintersectioncornerarea = getPolyCornerArea(opolys[x][y], cornerfacet[0], cornerfacet[1], cornerfacet[2])
                    totalunshiftedareas.append(abs(polyintersectioncornerarea))
                elif predfacettype == 'arc':
                    if predfacets[x][y][2] > 0:
                        #convex facet
                        polyintersectioncirclearea, _ = getCircleIntersectArea([predfacets[x][y][1][0], predfacets[x][y][1][1]], predfacets[x][y][2], opolys[x][y])
                    else:
                        #concave facet
                        polyintersectioncirclearea, _ = getCircleIntersectArea([predfacets[x][y][1][0], predfacets[x][y][1][1]], -predfacets[x][y][2], opolys[x][y])
                        polyintersectioncirclearea = getArea(opolys[x][y]) - polyintersectioncirclearea
                    totalunshiftedareas.append(abs(polyintersectioncirclearea))
                else: #curved corner:
                    curvedcorner = predfacets[x][y][-1]
                    polyintersectioncurvedcornerarea = getPolyCurvedCornerArea(opolys[x][y], curvedcorner[0], curvedcorner[1], curvedcorner[2], predfacets[x][y][3], predfacets[x][y][4])
                    totalunshiftedareas.append(abs(polyintersectioncurvedcornerarea))

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
                for testx in range(-checksize, checksize+1):
                    for testy in range(-checksize, checksize+1):
                        checkx = x-testx
                        checky = y-testy
                        if checkx >= 0 and checkx < len(opolys) and checky >= 0 and checky < len(opolys[0]):
                            testpoly = opolys[checkx][checky]
                            testbounds = [min(list(map(lambda x : x[0], testpoly))), min(list(map(lambda x : x[1], testpoly))), max(list(map(lambda x : x[0], testpoly))), max(list(map(lambda x : x[1], testpoly)))]
                            if not(testbounds[2] <= shiftbounds[0] or shiftbounds[2] <= testbounds[0] or testbounds[3] <= shiftbounds[1] or shiftbounds[3] <= testbounds[1]):
                                #bounding boxes intersect, could be nonzero intersection
                                polyintersections = getPolyIntersectArea(shiftpoly, testpoly)
                                if len(polyintersections) == 0:
                                    #No intersection
                                    continue
                                #For each overlap region
                                for polyintersection in polyintersections:
                                    if predfacettype == 'linear':
                                        shiftlinearfacet = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], predfacets[x][y][-1]))
                                        polyintersectionlineararea = getPolyLineArea(polyintersection, shiftlinearfacet[0], shiftlinearfacet[len(shiftlinearfacet)-1])
                                        nareas[checkx][checky].append([abs(polyintersectionlineararea), [x, y]])
                                    elif predfacettype == 'corner':
                                        shiftcornerfacet = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], predfacets[x][y][-1]))
                                        polyintersectioncornerarea = getPolyCornerArea(polyintersection, shiftcornerfacet[0], shiftcornerfacet[1], shiftcornerfacet[2])
                                        nareas[checkx][checky].append([abs(polyintersectioncornerarea), [x, y]])
                                    elif predfacettype == 'arc':
                                        if predfacets[x][y][2] > 0:
                                            polyintersectioncirclearea, _ = getCircleIntersectArea([predfacets[x][y][1][0]+velocity[0], predfacets[x][y][1][1]+velocity[1]], predfacets[x][y][2], polyintersection)
                                        else:
                                            polyintersectioncirclearea, _ = getCircleIntersectArea([predfacets[x][y][1][0]+velocity[0], predfacets[x][y][1][1]+velocity[1]], -predfacets[x][y][2], polyintersection)
                                            polyintersectioncirclearea = getArea(polyintersection) - polyintersectioncirclearea
                                        nareas[checkx][checky].append([abs(polyintersectioncirclearea), [x, y]])
                                    elif predfacettype == 'curvedcorner':
                                        curvedcorner = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], predfacets[x][y][-1]))
                                        polyintersectioncurvedcornerarea = getPolyCurvedCornerArea(polyintersection, curvedcorner[0], curvedcorner[1], curvedcorner[2], predfacets[x][y][3], predfacets[x][y][4])
                                        nareas[checkx][checky].append([abs(polyintersectioncurvedcornerarea), [x, y]])

            elif areas[x][y] == 1:
                shiftpoly = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], opolys[x][y]))
                shiftbounds = [min(list(map(lambda x : x[0], shiftpoly))), min(list(map(lambda x : x[1], shiftpoly))), max(list(map(lambda x : x[0], shiftpoly))), max(list(map(lambda x : x[1], shiftpoly)))]
                #bounds: minx, miny, maxx, maxy
                for testx in range(-checksize, checksize+1):
                    for testy in range(-checksize, checksize+1):
                        checkx = x-testx
                        checky = y-testy
                        if checkx >= 0 and checkx < len(opolys) and checky >= 0 and checky < len(opolys[0]):
                            testpoly = opolys[checkx][checky]
                            testbounds = [min(list(map(lambda x : x[0], testpoly))), min(list(map(lambda x : x[1], testpoly))), max(list(map(lambda x : x[0], testpoly))), max(list(map(lambda x : x[1], testpoly)))]
                            if not(testbounds[2] <= shiftbounds[0] or shiftbounds[2] <= testbounds[0] or testbounds[3] <= shiftbounds[1] or shiftbounds[3] <= testbounds[1]):
                                #bounding boxes intersect, could be nonzero intersection
                                polyintersections = getPolyIntersectArea(shiftpoly, testpoly)
                                if len(polyintersections) == 0:
                                    #No intersection
                                    continue
                                #For each overlap region
                                for polyintersection in polyintersections:
                                    nareas[checkx][checky].append([abs(getArea(polyintersection)), [x, y]])

            #Unpredicted partial fraction
            elif areas[x][y] > 0:
                print("Unpredicted partial fraction: {}, {}".format(x, y))

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