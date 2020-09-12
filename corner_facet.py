import numpy as np
import math
from geoms import getArea, lerp, getPolyLineArea, getPolyIntersectArea, getDistance
from circular_facet import getCircleIntersectArea, getCenter

#Returns area to left of corner facet1, corner, facet2
#TODO: more consistent way of computing polygon-corner intersection than extending corner? For example, test for polygon-line intersects.
def getPolyCornerArea(poly, facet1, corner, facet2):
    #Hyperparameters:
    extendcorner = 5 #Extend the sides of corners by this factor

    newfacet1 = lerp(facet1, corner, 1-extendcorner)
    newfacet2 = lerp(facet2, corner, 1-extendcorner)
    cornerpoly = [newfacet1, corner, newfacet2]
    cornerorientationarea = getArea(cornerpoly)
    if cornerorientationarea == 0:
        return 0
    elif cornerorientationarea > 0:
        #Convex vertex
        cornerintersect = getPolyIntersectArea(poly, cornerpoly)
        cornerareafraction = getPolyLineArea(poly, newfacet1, newfacet2)
        for cornerintersectpoly in cornerintersect:
            cornerareafraction += getArea(cornerintersectpoly)
    else:
        #Concave vertex
        reversecorner = cornerpoly.copy()
        reversecorner.reverse()
        cornerintersect = getPolyIntersectArea(poly, reversecorner)
        cornerareafraction = getPolyLineArea(poly, newfacet1, newfacet2)
        for cornerintersectpoly in cornerintersect:
            cornerareafraction -= getArea(cornerintersectpoly)
    return cornerareafraction

#Get region to left of line within polygon
def getPolyLineRegion(poly, l1, l2):
    assert not (l1[0] == l2[0] and l1[1] == l2[1])
    intersectRegion = []
    if l1[0] == l2[0]:
        for i in range(len(poly)):
            p1 = poly[i]
            p2 = poly[(i+1) % len(poly)]
            if (p1[0] <= l1[0] and l1[1] < l2[1]) or (p1[0] >= l1[0] and l1[1] > l2[1]):
                intersectRegion.append(p1)
            if (p1[0] < l1[0] and p2[0] > l1[0]) or (p1[0] > l1[0] and p2[0] < l1[0]):
                t = (l1[0] - p1[0])/(p2[0]-p1[0])
                pinter = lerp(p1, p2, t)
                intersectRegion.append(pinter)
    else:
        l = lambda x : l1[1] + (l2[1]-l1[1])*(x-l1[0])/(l2[0]-l1[0])
        for i in range(len(poly)):
            p1 = poly[i]
            p2 = poly[(i+1) % len(poly)]
            if (p1[1] >= l(p1[0]) and l1[0] < l2[0]) or (p1[1] <= l(p1[0]) and l1[0] > l2[0]):
                intersectRegion.append(p1)
            if (p1[1] > l(p1[0]) and p2[1] < l(p2[0])) or (p1[1] < l(p1[0]) and p2[1] > l(p2[0])):
                t = (p1[1]-l1[1]-(l2[1]-l1[1])*(p1[0]-l1[0])/(l2[0]-l1[0]))/((l2[1]-l1[1])*(p2[0]-p1[0])/(l2[0]-l1[0]) - (p2[1]-p1[1]))
                pinter = lerp(p1, p2, t)
                intersectRegion.append(pinter)
    return intersectRegion

#If radius1 is positive, center lies towards left of segment facet1 to corner; if negative, center on right
#If radius2 is positive, center lies towards left of segment corner to facet2; if negative, center on right
#If radius1 is None, treats segment facet1 to corner as a straight edge
#If radius2 is None, treats segment corner to facet2  as a straight edge
def getPolyCurvedCornerArea(poly, facet1, corner, facet2, radius1, radius2):
    if radius1 is not None:
        assert abs(radius1) >= getDistance(facet1, corner)/2
    if radius2 is not None:
        assert abs(radius2) >= getDistance(facet2, corner)/2
    area = getPolyCornerArea(poly, facet1, corner, facet2)*getArea(poly)
    if radius1 is not None:
        if radius1 < 0:
            radius1 *= -1
            polyfacet1region = getPolyLineRegion(poly, facet1, corner)
            facet1center = getCenter(corner, facet1, radius1)
            circlearea1, _ = getCircleIntersectArea(facet1center, radius1, polyfacet1region)
            circlearea1 *= -1
        else:
            polyfacet1region = getPolyLineRegion(poly, corner, facet1)
            facet1center = getCenter(facet1, corner, radius1)
            circlearea1, _ = getCircleIntersectArea(facet1center, radius1, polyfacet1region)
        area += circlearea1
    if radius2 is not None:
        if radius2 < 0:
            radius2 *= -1
            polyfacet2region = getPolyLineRegion(poly, corner, facet2)
            facet2center = getCenter(facet2, corner, radius2)
            circlearea2, _ = getCircleIntersectArea(facet2center, radius2, polyfacet2region)
            circlearea2 *= -1
        else:
            polyfacet2region = getPolyLineRegion(poly, facet2, corner)
            facet2center = getCenter(corner, facet2, radius2)
            circlearea2, _ = getCircleIntersectArea(facet2center, radius2, polyfacet2region)
        area += circlearea2
    return area

#Newton's method to optimize corner position on min least squares
#Matches area fractions in areas
#polys = list of grids; areas = corresponding list of area fractions
#epsilon = desired accuracy in area fractions
#convthreshold = alternate end condition, when changes in area fractions are all less than this threshold
#TODO: better way of testing convergence than current end condition with convthreshold?
def getCurvedCornerFacet(polys, areas, facet1, corner, facet2, radius1, radius2, epsilon, convthreshold, newtonfactor):
    #Hyperparameters:
    dtbase = 1e-8 #Timestep used to calculate numerical estimates of derivatives
    newtonFactor = 1 #Move this amount in Newton direction

    cornerx = corner[0]
    cornery = corner[1]
    curareas = list(map(lambda x : getPolyCurvedCornerArea(x, facet1, corner, facet2, radius1, radius2)/getArea(polys[x]), polys)) #area fractions
    while max(list(map(lambda x : abs(curareas[x] - areas[x]), range(len(polys))))) > epsilon:
        #Compute numerical derivatives
        curareasdx = list(map(lambda x : (getPolyCurvedCornerArea(x, facet1, [cornerx+dtbase/2, cornery], facet2, radius1, radius2)-getPolyCurvedCornerArea(x, facet1, [cornerx-dtbase/2, cornery], facet2, radius1, radius2))/dtbase, polys))
        curareasdy = list(map(lambda x : (getPolyCurvedCornerArea(x, facet1, [cornerx, cornery+dtbase/2], facet2, radius1, radius2)-getPolyCurvedCornerArea(x, facet1, [cornerx, cornery-dtbase/2], facet2, radius1, radius2))/dtbase, polys))
        #Jacobian: 2 by len(polys)
        jacobian = np.array([list(map(lambda x : curareasdx[x], range(len(polys)))), list(map(lambda x : curareasdy[x], range(len(polys))))])
        cornerchanges = np.matmul(np.matmul(np.linalg.inv(np.matmul(jacobian, np.transpose(jacobian))), jacobian), np.transpose(np.array(list(map(lambda x : areas[x]-curareas[x], range(len(polys)))))))
        cornerx += newtonFactor*cornerchanges[0]
        cornery += newtonFactor*cornerchanges[1]
        newcurareas = list(map(lambda x : getPolyCurvedCornerArea(x, facet1, [cornerx, cornery], facet2, radius1, radius2)/getArea(polys[x]), polys)) #area fractions
        if max(list(map(lambda x : abs(curareas[x]-newcurareas[x]), range(len(polys))))) < convthreshold:
            break
        curareas = newcurareas
        
    return [facet1, [cornerx, cornery], facet2]