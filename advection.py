import numpy as np
import math
import vtk
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from geoms import getArea, mergePolys, getPolyIntersectArea, getPolyLineArea, getPolyLineIntersects, lineIntersect
from linear_facet import getLinearFacet
from circular_facet import getCircleIntersectArea, getArcFacet
from mesh import makeCartesianGrid, makeQuadGrid, makeConcaveGrid

#opolys = grid of quads
#areas = grid of area fractions per quad
def merge(opolys, areas):
    #Merging
    dirs = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    mergedpolys = [] #coordinates of polys to be merged
    predmergedpolys = [[None] * len(opolys) for _ in range(len(opolys))]
    
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if areas[x][y] == 0 or areas[x][y] == 1 or predmergedpolys[x][y] is not None:
                continue
            #Compute number of neighbors
            gridsquare = [x, y]
            gridneighbors = []
            for direction in dirs:
                #Check directions in counterclockwise order
                testsquare = [gridsquare[0] + direction[0], gridsquare[1] + direction[1]]
                if testsquare[0] < len(opolys) and testsquare[0] >= 0 and testsquare[1] < len(opolys[0]) and testsquare[1] >= 0 and areas[testsquare[0]][testsquare[1]] < 1 and areas[testsquare[0]][testsquare[1]] > 0:
                    #This direction is a partial area neighbor
                    gridneighbors.append(testsquare)
            #Casework on number of neighbors
            #Case of single neighbor is eliminated in case of three neighbors
            if len(gridneighbors) == 2:
                #Path is clear
                predmergedpolys[x][y] = len(mergedpolys)
                mergedpolys.append([[gridsquare], gridneighbors])
            elif len(gridneighbors) >= 3:
                #Time to merge
                #invariant: mergeneighbors contains cells to be merged, whose neighbors may not have been explored yet; gridneighbors are the cells that need to be explored
                mergeneighbors = [gridsquare]
                for gridneighbor in gridneighbors:
                    mergeneighbors.append(gridneighbor)
                while len(gridneighbors) > 0:
                    testsquare = gridneighbors[0]
                    gridneighbors.pop(0)
                    testneighbors = []
                    for direction in dirs:
                        test2square = [testsquare[0] + direction[0], testsquare[1] + direction[1]]
                        if test2square[0] < len(opolys) and test2square[0] >= 0 and test2square[1] < len(opolys[0]) and test2square[1] >= 0 and areas[test2square[0]][test2square[1]] < 1 and areas[test2square[0]][test2square[1]] > 0:
                            testneighbors.append(test2square)
                    if len(testneighbors) >= 3:
                        #adjacent 3neighbor
                        for testneighbor in testneighbors:
                            if testneighbor not in mergeneighbors:
                                mergeneighbors.append(testneighbor)
                                gridneighbors.append(testneighbor)
                    elif len(testneighbors) == 2:
                        #test for case with two 3neighbors diagonal to each other
                        for testneighbor in testneighbors:
                            if testneighbor not in mergeneighbors:
                                testneighbor2count = 0
                                for direction in dirs:
                                    test2square = [testneighbor[0] + direction[0], testneighbor[1] + direction[1]]
                                    if test2square in mergeneighbors:
                                        testneighbor2count += 1
                                if testneighbor2count >= 2:
                                    mergeneighbors.append(testneighbor)
                                    gridneighbors.append(testneighbor)
                #mergeneighbors should consist of the cells to be merged + the two neighbors (not to be merged)
                mergeendpoints = []
                for mergeneighbor in mergeneighbors:
                    nummergeneighbors = 0
                    numallneighbors = 0
                    for direction in dirs:
                        test2square = [mergeneighbor[0] + direction[0], mergeneighbor[1] + direction[1]]
                        if test2square[0] < len(opolys) and test2square[0] >= 0 and test2square[1] < len(opolys[0]) and test2square[1] >= 0 and areas[test2square[0]][test2square[1]] < 1 and areas[test2square[0]][test2square[1]] > 0:
                            numallneighbors += 1
                            if test2square in mergeneighbors:
                                nummergeneighbors += 1
                    if nummergeneighbors < 2 and numallneighbors > 1:
                        #should be an endpoint
                        mergeendpoints.append(mergeneighbor)
                for mergeendpoint in mergeendpoints:
                    mergeneighbors.remove(mergeendpoint)
                for mergeneighbor in mergeneighbors:
                    predmergedpolys[mergeneighbor[0]][mergeneighbor[1]] = len(mergedpolys)
                mergedpolys.append([mergeneighbors, mergeendpoints])
    
    mergedcoords = [] #the coordinates of the fully merged poly
    mergedareafractions = [] #area fractions of fully merged poly
    for mergedpoly in mergedpolys:
        if len(mergedpoly[0]) == 1:
            mergedcoords.append(opolys[mergedpoly[0][0][0]][mergedpoly[0][0][1]])
            mergedareafractions.append(areas[mergedpoly[0][0][0]][mergedpoly[0][0][1]])
        else:
            boundary = opolys[mergedpoly[0][0][0]][mergedpoly[0][0][1]]
            for mergedpolyi in range(1, len(mergedpoly[0])):
                boundary = mergePolys(boundary, opolys[mergedpoly[0][mergedpolyi][0]][mergedpoly[0][mergedpolyi][1]])
            mergedcoords.append(boundary)
            boundaryarea = sum(list(map(lambda x : areas[x[0]][x[1]]*getArea(opolys[x[0]][x[1]]), mergedpoly[0])))/sum(list(map(lambda x : getArea(opolys[x[0]][x[1]]), mergedpoly[0])))
            mergedareafractions.append(boundaryarea)

    return predmergedpolys, mergedpolys, mergedcoords, mergedareafractions

#predmergedpolys, mergedpolys, mergedcoords, mergedareafractions: from merging function
#opolys, areas: original volume fractions, used to refine linear facets
#threshold: tolerance for errors in volume fractions
def makeFacets(predmergedpolys, mergedpolys, mergedcoords, mergedareafractions, opolys, areas, threshold):
    #facets
    predfacets = [[None] * len(predmergedpolys) for _ in range(len(predmergedpolys))]

    #orientations
    predorientations = [[None] * len(predmergedpolys) for _ in range(len(predmergedpolys))]
    mergedorientations = [None] * len(mergedpolys)
    
    #compute facet orientations
    for y in range(len(predmergedpolys[0])):
        for x in range(len(predmergedpolys)):
            if predmergedpolys[x][y] is not None and mergedorientations[predmergedpolys[x][y]] == None:
                #finding path, starting with bottom left most point
                #orientation of path: outer boundary = True, inner boundary = False
                if (y > 0 and areas[x][y-1] == 1) or (x > 0 and areas[x-1][y] == 1):
                    curOrientation = False
                else:
                    curOrientation = True
                path = []
                curmergedpolyindex = predmergedpolys[x][y]
                curmergedpoly = mergedpolys[curmergedpolyindex]
                if curmergedpoly[1][0][0] < curmergedpoly[1][1][0]:
                    nextmergedpolyindex = predmergedpolys[curmergedpoly[1][1][0]][curmergedpoly[1][1][1]]
                elif curmergedpoly[1][0][0] > curmergedpoly[1][1][0]:
                    nextmergedpolyindex = predmergedpolys[curmergedpoly[1][0][0]][curmergedpoly[1][0][1]]
                else:
                    if curmergedpoly[1][0][1] < curmergedpoly[1][1][1]:
                        nextmergedpolyindex = predmergedpolys[curmergedpoly[1][0][0]][curmergedpoly[1][0][1]]
                    else:
                        nextmergedpolyindex = predmergedpolys[curmergedpoly[1][1][0]][curmergedpoly[1][1][1]]
                path.append(curmergedpolyindex)
                mergedorientations[curmergedpolyindex] = curOrientation
                for point in curmergedpoly[0]:
                    predorientations[point[0]][point[1]] = curOrientation
                prevmergedpolyindex = curmergedpolyindex
                while nextmergedpolyindex != path[0]:
                    curmergedpolyindex = nextmergedpolyindex
                    curmergedpoly = mergedpolys[curmergedpolyindex]
                    if curmergedpoly[1][0] in mergedpolys[prevmergedpolyindex][0]:
                        nextmergedpolyindex = predmergedpolys[curmergedpoly[1][1][0]][curmergedpoly[1][1][1]]
                    else:
                        assert curmergedpoly[1][1] in mergedpolys[prevmergedpolyindex][0]
                        nextmergedpolyindex = predmergedpolys[curmergedpoly[1][0][0]][curmergedpoly[1][0][1]]
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
                    #grid1, 2, 3 are indices of mergedpolys
                    grid1 = path[pathelement-1]
                    grid2 = path[pathelement]
                    grid3 = path[pathelement+1]
                    try:
                        if mergedorientations[grid2]:
                            prevpolygon = mergedcoords[grid1]
                            prevpolygonarea = mergedareafractions[grid1]*getArea(mergedcoords[grid1])
                            nextpolygon = mergedcoords[grid3]
                            nextpolygonarea = mergedareafractions[grid3]*getArea(mergedcoords[grid3])
                        else:
                            prevpolygon = mergedcoords[grid1]
                            prevpolygonarea = (1-mergedareafractions[grid1])*getArea(mergedcoords[grid1])
                            nextpolygon = mergedcoords[grid3]
                            nextpolygonarea = (1-mergedareafractions[grid3])*getArea(mergedcoords[grid3])
                        facetline1, facetline2 = getLinearFacet(prevpolygon, nextpolygon, prevpolygonarea, nextpolygonarea, threshold)
                        if abs(mergedareafractions[grid2] - getPolyLineArea(mergedcoords[grid2], facetline1, facetline2)/getArea(mergedcoords[grid2])) > threshold:
                            #maybe excess merging reduced resolution: try immediate neighbors
                            if mergedpolys[grid2][1][0] in mergedpolys[grid1][0]:
                                prevneighbor = mergedpolys[grid2][1][0]
                                nextneighbor = mergedpolys[grid2][1][1]
                            else:
                                assert mergedpolys[grid2][1][1] in mergedpolys[grid1][0]
                                prevneighbor = mergedpolys[grid2][1][1]
                                nextneighbor = mergedpolys[grid2][1][0]
                                
                            if mergedorientations[grid2]:
                                prevpolygon = opolys[prevneighbor[0]][prevneighbor[1]]
                                prevpolygonarea = areas[prevneighbor[0]][prevneighbor[1]]*getArea(opolys[prevneighbor[0]][prevneighbor[1]])
                                nextpolygon = opolys[nextneighbor[0]][nextneighbor[1]]
                                nextpolygonarea = areas[nextneighbor[0]][nextneighbor[1]]*getArea(opolys[nextneighbor[0]][nextneighbor[1]])
                            else:
                                prevpolygon = opolys[prevneighbor[0]][prevneighbor[1]]
                                prevpolygonarea = (1-areas[prevneighbor[0]][prevneighbor[1]])*getArea(opolys[prevneighbor[0]][prevneighbor[1]])
                                nextpolygon = opolys[nextneighbor[0]][nextneighbor[1]]
                                nextpolygonarea = (1-areas[nextneighbor[0]][nextneighbor[1]])*getArea(opolys[nextneighbor[0]][nextneighbor[1]])
                            facetline1, facetline2 = getLinearFacet(prevpolygon, nextpolygon, prevpolygonarea, nextpolygonarea, threshold)
                        if abs(mergedareafractions[grid2] - getPolyLineArea(mergedcoords[grid2], facetline1, facetline2)/getArea(mergedcoords[grid2])) < threshold*10:
                            intersects = getPolyLineIntersects(mergedcoords[grid2], facetline1, facetline2)
                            for linearfacetsquares in mergedpolys[grid2][0]:
                                predfacets[linearfacetsquares[0]][linearfacetsquares[1]] = ['linear', intersects]
                            facetfitted[pathelement-1] = intersects
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
                        curlinearfacet = facetfitted[(pathelement) % len(facetfitted)]
                for pathelement in range(1, 2*len(facetfitted)-1):
                    if facetfitted[(len(facetfitted)-1-pathelement) % len(facetfitted)] is None:
                        rights[(len(facetfitted)-1-pathelement) % len(facetfitted)] = curlinearfacet
                    else:
                        curlinearfacet = facetfitted[(len(facetfitted)-1-pathelement) % len(facetfitted)]
                        
                #Try to make corners
                #here pathelement is 1 less than pathelement from linear facet fitting
                for pathelement in range(len(facetfitted)):
                    if facetfitted[pathelement] is None and lefts[pathelement] is not None and rights[pathelement] is not None and lefts[pathelement] != rights[pathelement]:
                        prevFacet = lefts[pathelement]
                        nextFacet = rights[pathelement]
                        corner = lineIntersect(prevFacet[0], prevFacet[1], nextFacet[0], nextFacet[1])
                        corner = [prevFacet[0], corner, nextFacet[1]]
                            #Convex vertex
                        if getArea(corner) > 0:
                            cornerintersect = getPolyIntersectArea(mergedcoords[path[pathelement+1]], corner)
                            cornerareafraction = getPolyLineArea(mergedcoords[path[pathelement+1]], prevFacet[0], nextFacet[1])
                            for cornerintersectpoly in cornerintersect:
                                cornerareafraction += getArea(cornerintersectpoly)
                            cornerareafraction /= getArea(mergedcoords[path[pathelement+1]])
                        #Concave vertex
                        else:
                            reversecorner = corner.copy()
                            reversecorner.reverse()
                            cornerintersect = getPolyIntersectArea(mergedcoords[path[pathelement+1]], reversecorner)
                            cornerareafraction = getPolyLineArea(mergedcoords[path[pathelement+1]], prevFacet[0], nextFacet[1])
                            for cornerintersectpoly in cornerintersect:
                                cornerareafraction -= getArea(cornerintersectpoly)
                            cornerareafraction /= getArea(mergedcoords[path[pathelement+1]])
                        #Invert if inverted
                        if not(mergedorientations[path[pathelement+1]]):
                            cornerareafraction = 1-cornerareafraction
                        #Make corner if close enough
                        if abs(cornerareafraction - mergedareafractions[path[pathelement+1]]) < threshold*1e2:
                            #print("Corner facet: {}".format(path[pathelement+1]))
                            for cornerfacetsquares in mergedpolys[path[pathelement+1]][0]:
                                predfacets[cornerfacetsquares[0]][cornerfacetsquares[1]] = ['corner', corner]
                            facetfitted[pathelement] = corner
                        else:
                            continue
                                
                #Try to make arc facets for the remaining ones
                for pathelement in range(1, len(path)-1):
                    if facetfitted[pathelement-1] is None:
                        #grid1, 2, 3 are indices of mergedpolys
                        grid1 = path[pathelement-1]
                        grid2 = path[pathelement]
                        grid3 = path[pathelement+1]
                        try:
                            if mergedorientations[grid2]:
                                prevpolygon = mergedcoords[grid1]
                                prevpolygonarea = mergedareafractions[grid1]
                                curpolygonarea = mergedareafractions[grid2]
                                nextpolygon = mergedcoords[grid3]
                                nextpolygonarea = mergedareafractions[grid3]
                            else:
                                prevpolygon = mergedcoords[grid1]
                                prevpolygonarea = 1-mergedareafractions[grid1]
                                curpolygonarea = 1-mergedareafractions[grid2]
                                nextpolygon = mergedcoords[grid3]
                                nextpolygonarea = 1-mergedareafractions[grid3]
                            arccenter, arcradius, arcintersects = getArcFacet(prevpolygon, mergedcoords[grid2], nextpolygon, prevpolygonarea, curpolygonarea, nextpolygonarea, threshold)
                            for arcfacetsquares in mergedpolys[grid2][0]:
                                predfacets[arcfacetsquares[0]][arcfacetsquares[1]] = ['arc', arccenter, arcradius, arcintersects]
                            facetfitted[pathelement-1] = [arccenter, arcradius, arcintersects]
                        except:
                            continue

    return predfacets, predorientations

def advectFacets(opolys, areas, predfacets, predorientations, velocity, threshold):
    #areas after new timestep
    nareas = [[0] * (len(opolys)) for _ in range(len(opolys))]
    
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            if predfacets[x][y] is not None:
                predfacettype = predfacets[x][y][0]
                shiftpoly = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], opolys[x][y]))
                shiftbounds = [min(list(map(lambda x : x[0], shiftpoly))), min(list(map(lambda x : x[1], shiftpoly))), max(list(map(lambda x : x[0], shiftpoly))), max(list(map(lambda x : x[1], shiftpoly)))]
                #bounds: minx, miny, maxx, maxy
                #For each neighbor of the cell
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
                                                nareas[x-testx][y-testy] += abs(polyintersectionlineararea)
                                            else:
                                                nareas[x-testx][y-testy] += abs(getArea(polyintersection))-abs(polyintersectionlineararea)
                                        elif predfacettype == 'corner':
                                            shiftcornerfacet = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], predfacets[x][y][1]))
                                            if getArea(shiftcornerfacet) > 0:
                                                cornerintersections = getPolyIntersectArea(shiftcornerfacet, polyintersection)
                                                polyintersectioncornerarea = abs(getPolyLineArea(polyintersection, shiftcornerfacet[0], shiftcornerfacet[len(shiftcornerfacet)-1]))
                                                for cornerintersection in cornerintersections:
                                                    polyintersectioncornerarea += abs(getArea(cornerintersection))
                                            else:
                                                shiftcornerfacet.reverse()
                                                cornerintersections = getPolyIntersectArea(shiftcornerfacet, polyintersection)
                                                polyintersectioncornerarea = abs(getPolyLineArea(polyintersection, shiftcornerfacet[len(shiftcornerfacet)-1], shiftcornerfacet[0]))
                                                for cornerintersection in cornerintersections:
                                                    polyintersectioncornerarea -= abs(getArea(cornerintersection))
                                            if predorientations[x][y]:
                                                nareas[x-testx][y-testy] += abs(polyintersectioncornerarea)
                                            else:
                                                nareas[x-testx][y-testy] += abs(getArea(polyintersection))-abs(polyintersectioncornerarea)
                                        else: #circular
                                            polyintersectioncirclearea, _ = getCircleIntersectArea([predfacets[x][y][1][0]+velocity[0], predfacets[x][y][1][1]+velocity[1]], predfacets[x][y][2], polyintersection)
                                            if predorientations[x][y]:
                                                nareas[x-testx][y-testy] += abs(polyintersectioncirclearea)
                                            else:
                                                nareas[x-testx][y-testy] += abs(getArea(polyintersection))-abs(polyintersectioncirclearea)
                                except:
                                    print("Error in computing {} intersection of mesh cells: ({}, {}), ({}, {})".format(predfacettype, x, y, x-testx, y-testy))
            elif areas[x][y] == 1:
                shiftpoly = list(map(lambda x : [x[0]+velocity[0], x[1]+velocity[1]], opolys[x][y]))
                shiftbounds = [min(list(map(lambda x : x[0], shiftpoly))), min(list(map(lambda x : x[1], shiftpoly))), max(list(map(lambda x : x[0], shiftpoly))), max(list(map(lambda x : x[1], shiftpoly)))]
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
                                        nareas[x-testx][y-testy] += abs(getArea(polyintersection))
                                except:
                                    print("Error in computing intersection of full mesh cells: ({}, {}), ({}, {})".format(x, y, x-testx, y-testy))

            #Unpredicted partial fraction
            elif areas[x][y] > 0:
                print("Unpredicted partial fraction: {}, {}".format(x, y))
                                    
    #convert areas into area fractions, fix errors with area intersection computation
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            newarea = nareas[x][y]/getArea(opolys[x][y])
            if abs(1 - newarea) < threshold:
                nareas[x][y] = 1
            elif abs(newarea) < threshold:
                nareas[x][y] = 0
            else:
                nareas[x][y] = newarea

    return nareas

def advection(gridSize, timesteps):
    
    #starting mesh
    wiggle = 0.50
    opoints = makeCartesianGrid(gridSize)
    #opoints = makeQuadGrid(gridSize, wiggle)
    #opoints = makeConcaveGrid(gridSize, wiggle)
    
    #mesh velocity per timestep
    vsize = 0.095
    velocity = [vsize, vsize]
    threshold = 1e-10
    
    #variables per time step
    opolys = [[0] * (len(opoints)-1) for _ in range(len(opoints)-1)]
                
    #vtk
    opointindices = [[-1] * (len(opolys)+1) for _ in range(len(opolys)+1)]
    oindexcounter = 0
    ovtkpoints = vtk.vtkPoints()
    ougrid = vtk.vtkUnstructuredGrid()
    ougrid.Allocate((len(opolys)+1)**2)

    #matplotlib
    opatches = []
    opatchareas = []
    opatchareapartials = []
    
    #compute initial area setup
    areas = [[0] * len(opolys) for _ in range(len(opolys))]
    intersects = [[0] * len(opolys) for _ in range(len(opolys))]

    for x in range(len(opoints)-1):
        for y in range(len(opoints)-1):
            opoly = [opoints[x][y], opoints[x+1][y], opoints[x+1][y+1], opoints[x][y+1]]
            neighbors = [[0, 0], [1, 0], [1, 1], [0, 1]]
            for neighbor in neighbors:
                if opointindices[x+neighbor[0]][y+neighbor[1]] == -1:
                    opointindices[x+neighbor[0]][y+neighbor[1]] = oindexcounter
                    ovtkpoints.InsertPoint(oindexcounter, [opoints[x+neighbor[0]][y+neighbor[1]][0], opoints[x+neighbor[0]][y+neighbor[1]][1], 0])
                    oindexcounter += 1
            ougrid.InsertNextCell(vtk.VTK_QUAD, 4, [opointindices[x][y], opointindices[x+1][y], opointindices[x+1][y+1], opointindices[x][y+1]])
            
            #Initial area setting
            #"""
            radius = 6
            xpoints = [[radius/math.sqrt(2), 0], [2*radius/math.sqrt(2), radius/math.sqrt(2)], [3*radius/math.sqrt(2), 0], [4*radius/math.sqrt(2), radius/math.sqrt(2)], [3*radius/math.sqrt(2), 2*radius/math.sqrt(2)], [4*radius/math.sqrt(2), 3*radius/math.sqrt(2)], [3*radius/math.sqrt(2), 4*radius/math.sqrt(2)], [2*radius/math.sqrt(2), 3*radius/math.sqrt(2)], [radius/math.sqrt(2), 4*radius/math.sqrt(2)], [0, 3*radius/math.sqrt(2)], [radius/math.sqrt(2), 2*radius/math.sqrt(2)], [0, radius/math.sqrt(2)]]
            xpoints = list(map(lambda x : [x[0]+1.05, x[1]+1.05], xpoints))
            xpolyintersects = getPolyIntersectArea(xpoints, opoly)
            for xpolyintersect in xpolyintersects:
                areas[x][y] += abs(getArea(xpolyintersect))
            #"""
            """
            radius = 6
            radiussmall = 3
            center = [1.5*radius+0.025, 1.5*radius+0.025]
            area, intersect = getCircleIntersectArea(center, radius, opoly)
            areas[x][y] = area
            area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
            areas[x][y] -= area
            """
            areas[x][y] /= getArea(opoly)
            if abs(1-areas[x][y]) < threshold:
                areas[x][y] = 1
            elif abs(areas[x][y]) < threshold:
                areas[x][y] = 0
            
            #matplotlib: match partial areas and polygons for plotting
            opolys[x][y] = opoly
            opatch = Polygon(np.array(opoly), True)
            opatches.append(opatch)
            opatchareas.append(areas[x][y])
            opatchareapartials.append(math.ceil(areas[x][y] - math.floor(areas[x][y])))
                
    #plot in matplotlib
    op = PatchCollection(opatches, cmap='jet')
    opatchareas = np.array(opatchareas)
    op.set_array(opatchareas)
    fig, ax = plt.subplots()
    ax.set_xlim(-1, len(opolys)+1)
    ax.set_ylim(-1, len(opolys)+1)
    ax.add_collection(op)
    plt.savefig("advection_original.png", dpi=199)
    plt.clf()
    
    op = PatchCollection(opatches, cmap='jet')
    opatchareas = np.array(opatchareapartials)
    op.set_array(opatchareas)
    fig, ax = plt.subplots()
    ax.set_xlim(-1, len(opolys)+1)
    ax.set_ylim(-1, len(opolys)+1)
    ax.add_collection(op)
    plt.savefig("advection_original_partials.png", dpi=199)
    plt.clf()
    
    #plot in vtk
    plotpolys = []
    for x in range(len(opolys)):
        for y in range(len(opolys)):
            plotpolys.append(opolys[x][y])
    
    #for histogram
    prederrors = []
    
    #main advection loop
    for timestep in range(timesteps):
        print("Timestep: {}".format(timestep))
        
        predmergedpolys, mergedpolys, mergedcoords, mergedareafractions = merge(opolys, areas)
                
        predfacets, predorientations = makeFacets(predmergedpolys, mergedpolys, mergedcoords, mergedareafractions, opolys, areas, threshold)

        #Plot facets in vtk
        
        print("Computing new areas")
        nareas = advectFacets(opolys, areas, predfacets, predorientations, velocity, threshold)
        areas = nareas

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
        ax.set_xlim(-1, len(opolys)+1)
        ax.set_ylim(-1, len(opolys)+1)
        ax.add_collection(newp)
        plt.savefig("advection_pred_timestep_{}.png".format(timestep), dpi=199)
        plt.clf()

        newp = PatchCollection(opatches, cmap='jet')
        npatchareas = []
        for x in range(len(opolys)):
            for y in range(len(opolys)):
                npatchareas.append(math.ceil(areas[x][y] - math.floor(areas[x][y])))
        npatchareas = np.array(npatchareas)
        newp.set_array(npatchareas)
        fig, ax = plt.subplots()
        ax.set_xlim(-1, len(opolys)+1)
        ax.set_ylim(-1, len(opolys)+1)
        ax.add_collection(newp)
        plt.savefig("advection_predpartials_timestep_{}.png".format(timestep), dpi=199)
        plt.clf()
        
        plt.close("all")

np.random.seed(17)
gridSize = 50
timesteps = 10
advection(gridSize, timesteps)