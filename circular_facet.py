import math
import numpy as np
from geoms import getArea, getDistance, lerp, getPolyLineArea, getPolyLineIntersects, getCentroid, pointInPoly
from linear_facet import getLinearFacet

#Assume points are actually on circle
def getArcArea(p1, p2, center, radius):
    #Compute area
    distance = getDistance(p1, p2)
    theta = math.asin(distance/(2*radius))
    area = (theta * radius**2) - (distance * math.cos(theta)*radius / 2)
    #Check orientation of segment p1p2 compared to center
    if (p1[0]-center[0])*(p2[1]-p1[1]) - (p1[1]-center[1])*(p2[0]-p1[0]) >= 0:
        #Angle from p1 to p2 moves counterclockwise: return the small area
        return area
    #Return the large area (mod pi*radius^2)
    return -area

#Get center to the left of line l1 to l2
def getCenter(l1, l2, radius):
    assert getDistance(l1, l2) <= 2*radius, "Error in getCenter arguments: radius too small"
    center = lerp(l1, l2, 0.5)
    c = math.sqrt(radius**2 - (getDistance(l1, l2)**2)/4)/getDistance(l1, l2)
    return [center[0]+(c*(l1[1]-l2[1])), center[1]+(c*(l2[0]-l1[0]))]

#Get circle-circle intersects
def getCircleCircleIntersects(center1, center2, radius1, radius2):
    d = getDistance(center1, center2)
    if abs(radius1-radius2) > d or radius1+radius2 < d or d == 0:
        return []
    else:
        chorddistance = (abs(radius1**2 - radius2**2)+d**2)/(2*d)
        halfchordlength = math.sqrt(radius1**2 - chorddistance**2)
        chordmidpoint = lerp(center1, center2, chorddistance/d)
        return [[chordmidpoint[0]+halfchordlength*(center2[1]-center1[1])/d, chordmidpoint[1]-halfchordlength*(center2[0]-center1[0])/d], [chordmidpoint[0]-halfchordlength*(center2[1]-center1[1])/d, chordmidpoint[1]+halfchordlength*(center2[0]-center1[0])/d]]

#Line l1 to l2
def getCircleLineIntersects(l1, l2, center, radius):
    intersects = []
    a = (l2[0]-l1[0])**2 + (l2[1]-l1[1])**2
    b = 2*((l1[0]-center[0])*(l2[0]-l1[0]) + (l1[1]-center[1])*(l2[1]-l1[1]))
    c = (l1[0]-center[0])**2 + (l1[1]-center[1])**2 - radius**2
    #a is never 0 because p1 != p2
    disc = b**2 - 4*a*c
    if disc > 0:
        x1 = (-b - math.sqrt(disc))/(2*a)
        if x1 <= 1 and x1 > 0:
            inter1 = lerp(l1, l2, x1)
            intersects.append(inter1)
        x2 = (-b + math.sqrt(disc))/(2*a)
        if x2 < 1 and x2 >= 0:
            inter2 = lerp(l1, l2, x2)
            intersects.append(inter2)
    return intersects

#Line l1 to l2
def getCircleLineIntersects2(l1, l2, center, radius):
    intersects = []
    a = (l2[0]-l1[0])**2 + (l2[1]-l1[1])**2
    b = 2*((l1[0]-center[0])*(l2[0]-l1[0]) + (l1[1]-center[1])*(l2[1]-l1[1]))
    c = (l1[0]-center[0])**2 + (l1[1]-center[1])**2 - radius**2
    #a is never 0 because p1 != p2
    disc = b**2 - 4*a*c
    if disc > 0:
        x1 = (-b - math.sqrt(disc))/(2*a)
        inter1 = lerp(l1, l2, x1)
        intersects.append(inter1)
        x2 = (-b + math.sqrt(disc))/(2*a)
        inter2 = lerp(l1, l2, x2)
        intersects.append(inter2)
    return intersects

def getCircleIntersectArea(center, radius, poly):
    #Hyperparameter
    adjustcorneramount = 1e-14
    notmod = True
    while notmod:
        startAt = 1
        intersectpoints = []
        arcpoints = []
        for i in range(len(poly)):
            curpoint = poly[i]
            nextpoint = poly[(i+1) % len(poly)]
            curin = getDistance(curpoint, center) <= radius
            nextin = getDistance(nextpoint, center) <= radius
            if curin:
                intersectpoints.append(curpoint)
            lineintersects = getCircleLineIntersects(curpoint, nextpoint, center, radius)
            for intersect in lineintersects:
                intersectpoints.append(intersect)
                if len(arcpoints) == 0 and curin and not(nextin):
                    startAt = 0
                arcpoints.append(intersect)
        #If not 0 mod 2, circle intersects a corner, perturb poly and rerun
        if len(arcpoints) % 2 == 1:
            poly = list(map(lambda x : [x[0]+adjustcorneramount, x[1]+adjustcorneramount], poly))
        else:
            notmod = False

    area = 0
    for i in range(0, len(arcpoints), 2):
        area += getArcArea(arcpoints[startAt+i], arcpoints[(startAt+i+1) % len(arcpoints)], center, radius)
    area += getArea(intersectpoints)
    
    if len(arcpoints) == 0 and area == 0:
        if pointInPoly(center, poly):
            #circle lies entirely inside poly
            return radius*radius*math.pi, []
        #circle lies entirely outside poly
        return 0, []
    
    return area % (radius*radius*math.pi), arcpoints

#Newton's method to match a1, a2, a3s
#Matches area fractions a1, a2, a3
def getArcFacet(poly1, poly2, poly3, a1, a2, a3, epsilon):
    #Basic sanity tests:
    assert a1 >= 0 and a1 <= 1 and a2 >= 0 and a2 <= 1 and a3 >= 0 and a3 <= 1, print("Given areas for circular facet are not valid")

    #Hyperparameters:
    scaleEpsilon = 0.99 #Adjust threshold by this amount per failed timestep
    dtbase = 1e-8 #Timestep used to calculate numerical estimates of derivatives
    fixLinearFacetOrientation = epsilon #Threshold used to determine orientation of wrong linear facets
    fixLinearFacet = 1.1 #Amount to adjust linear facet guess by to fix
    maxTimestep = 500 #Amount of timesteps allowed before we declare failure
    
    #Rotate so that x-axis is linear facet
    l1, l2 = getLinearFacet(poly1, poly3, a1, a3, epsilon/10)

    #Center to left of line l1 to l2
    if getPolyLineArea(poly2, l1, l2) > a2*getArea(poly2):
        retcenter, retradius, retintersects = getArcFacet(poly3, poly2, poly1, 1-a3, 1-a2, 1-a1, epsilon)
        if retcenter is not None:
            retintersects.reverse()
            retradius *= -1
            return retcenter, retradius, retintersects #TODO: radius should always be positive, negative radius means facet is concave, area fraction k -> 1-k
        else:
            return None, None, None

    #Convert area fractions to areas
    poly1area = getArea(poly1)
    poly2area = getArea(poly2)
    poly3area = getArea(poly3)
    a1 *= poly1area
    a2 *= poly2area
    a3 *= poly3area

    poly1intersects = getPolyLineIntersects(poly1, l1, l2)
    poly2intersects = getPolyLineIntersects(poly2, l1, l2)
    poly3intersects = getPolyLineIntersects(poly3, l1, l2)

    if len(poly2intersects) > 0 and (getDistance(poly2intersects[0], poly1intersects[-1]) > fixLinearFacetOrientation and getDistance(poly2intersects[-1], poly3intersects[0]) > fixLinearFacetOrientation):
        #Not a good linear facet
        if getDistance(poly2intersects[-1], poly1intersects[0]) < fixLinearFacetOrientation:
            #Linear facet intersects poly2, poly1, poly3
            invertpoint = poly2intersects[-1]
            for p in poly2:
                if p in poly1:
                    invertpoint2 = p
                    break
            nl0 = l1[0] + fixLinearFacet*((invertpoint[0]-l1[0]) - ((invertpoint[0]-l1[0])*(invertpoint2[0]-invertpoint[0]) + (invertpoint[1]-l1[1])*(invertpoint2[1]-invertpoint[1]))/getDistance(invertpoint, invertpoint2)**2*(invertpoint2[0]-invertpoint[0]))
            nl1 = l1[1] + fixLinearFacet*((invertpoint[1]-l1[1]) - ((invertpoint[0]-l1[0])*(invertpoint2[0]-invertpoint[0]) + (invertpoint[1]-l1[1])*(invertpoint2[1]-invertpoint[1]))/getDistance(invertpoint, invertpoint2)**2*(invertpoint2[1]-invertpoint[1]))
            l1 = [nl0, nl1]
            nl0 = l2[0] + fixLinearFacet*((invertpoint[0]-l2[0]) - ((invertpoint[0]-l2[0])*(invertpoint2[0]-invertpoint[0]) + (invertpoint[1]-l2[1])*(invertpoint2[1]-invertpoint[1]))/getDistance(invertpoint, invertpoint2)**2*(invertpoint2[0]-invertpoint[0]))
            nl1 = l2[1] + fixLinearFacet*((invertpoint[1]-l2[1]) - ((invertpoint[0]-l2[0])*(invertpoint2[0]-invertpoint[0]) + (invertpoint[1]-l2[1])*(invertpoint2[1]-invertpoint[1]))/getDistance(invertpoint, invertpoint2)**2*(invertpoint2[1]-invertpoint[1]))
            l2 = [nl0, nl1]
        elif getDistance(poly2intersects[0], poly3intersects[-1]) < fixLinearFacetOrientation:
            #Linear facet intersects poly1, poly3, poly2
            invertpoint = poly2intersects[0]
            for p in poly2:
                if p in poly3:
                    invertpoint2 = p
                    break
            nl0 = l1[0] + fixLinearFacet*((invertpoint[0]-l1[0]) - ((invertpoint[0]-l1[0])*(invertpoint2[0]-invertpoint[0]) + (invertpoint[1]-l1[1])*(invertpoint2[1]-invertpoint[1]))/getDistance(invertpoint, invertpoint2)**2*(invertpoint2[0]-invertpoint[0]))
            nl1 = l1[1] + fixLinearFacet*((invertpoint[1]-l1[1]) - ((invertpoint[0]-l1[0])*(invertpoint2[0]-invertpoint[0]) + (invertpoint[1]-l1[1])*(invertpoint2[1]-invertpoint[1]))/getDistance(invertpoint, invertpoint2)**2*(invertpoint2[1]-invertpoint[1]))
            l1 = [nl0, nl1]
            nl0 = l2[0] + fixLinearFacet*((invertpoint[0]-l2[0]) - ((invertpoint[0]-l2[0])*(invertpoint2[0]-invertpoint[0]) + (invertpoint[1]-l2[1])*(invertpoint2[1]-invertpoint[1]))/getDistance(invertpoint, invertpoint2)**2*(invertpoint2[0]-invertpoint[0]))
            nl1 = l2[1] + fixLinearFacet*((invertpoint[1]-l2[1]) - ((invertpoint[0]-l2[0])*(invertpoint2[0]-invertpoint[0]) + (invertpoint[1]-l2[1])*(invertpoint2[1]-invertpoint[1]))/getDistance(invertpoint, invertpoint2)**2*(invertpoint2[1]-invertpoint[1]))
            l2 = [nl0, nl1]
        else:
            #Unknown configuration
            print("Unknown arc facet configuration")

    rot = lambda p : [((l2[0]-l1[0])*(p[0]-l1[0]) + (l2[1]-l1[1])*(p[1]-l1[1]))/getDistance(l1, l2), ((l1[1]-l2[1])*(p[0]-l1[0]) + (l2[0]-l1[0])*(p[1]-l1[1]))/getDistance(l1, l2)]
    unrot = lambda p : [((l2[0]-l1[0])*p[0] - (l2[1]-l1[1])*p[1])/getDistance(l1, l2) + l1[0], (-(l1[1]-l2[1])*p[0] + (l2[0]-l1[0])*p[1])/getDistance(l1, l2) + l1[1]]
    rpoly1 = list(map(rot, poly1))
    rpoly2 = list(map(rot, poly2))
    rpoly3 = list(map(rot, poly3))
    
    #Working in rotated frame
    #Find where x-axis intersects rpoly1, rpoly3
    intersects1 = []
    for i in range(len(rpoly1)):
        p1 = rpoly1[i]
        p2 = rpoly1[(i+1) % len(rpoly1)]
        if ((p1[1] > 0) != (p2[1] > 0)):
            intersects1.append(p1[0] + (p2[0]-p1[0])*p1[1]/(p1[1]-p2[1]))
    #towards r1down = cura1 decreases
    r1down = max(intersects1)
    r1up = min(intersects1)
    intersects3 = []
    for i in range(len(rpoly3)):
        p1 = rpoly3[i]
        p2 = rpoly3[(i+1) % len(rpoly3)]
        if ((p1[1] > 0) != (p2[1] > 0)):
            intersects3.append(p1[0] + (p2[0]-p1[0])*p1[1]/(p1[1]-p2[1]))
    #towards r3down = cura3 decreases
    r3down = min(intersects3)
    r3up = max(intersects3)

    #Idea: degrees of freedom are a, r1, r3
    #r1 lies between r1down and r1up, r3 lies between r3down and r3up
    #Adjust a to match a2, then adjust r1 to match a1, then adjust r3 to match a3, then repeat until convergence
    cura1 = float("inf")
    cura2 = float("inf")
    cura3 = float("inf")
    #higher t1, t3 = lower cura1, cura3 respectively
    t1 = 0.5
    t3 = 0.5
    converged = False
    r1 = r1down*t1 + r1up*(1-t1)
    r3 = r3down*t3 + r3up*(1-t3)
    radius = abs(r3-r1)
    center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
    
    numcycles = 0
    
    while (not(converged) and numcycles < maxTimestep):

        #adjust radius to match a2
        radius = abs(r3-r1)
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura2, _ = getCircleIntersectArea(center, radius, rpoly2)
        rgap = radius
        doConverge = False
        #higher a = higher cura2
        largeenough = False
        while abs(cura2 - a2)/poly2area > epsilon*(scaleEpsilon**numcycles):
            if not(largeenough):
                largeenough2 = True
                for rpoly2vertex in rpoly2:
                    if rpoly2vertex[1] > 0 and getDistance(rpoly2vertex, center) > radius:
                        largeenough2 = False
                if largeenough2:
                    largeenough = True
            if cura2 > a2:
                radius += rgap
                if not(doConverge):
                    rgap *= 2
                else:
                    rgap /= 2
            else:
                if not(largeenough):
                    radius += rgap
                    rgap *= 2
                elif not(doConverge):
                    rgap /= 4
                    radius -= rgap
                    doConverge = True
                    rgap /= 2
                else:
                    radius -= rgap
                    rgap /= 2
            center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
            cura2, _ = getCircleIntersectArea(center, radius, rpoly2)

        #adjust r1 to match a1
        dt = dtbase
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura1, _ = getCircleIntersectArea(center, radius, rpoly1)
        
        while abs(cura1 - a1)/poly1area > epsilon*(scaleEpsilon**numcycles):
            cura1plusdt, _ = getCircleIntersectArea([(r1+dt+r3)/2, math.sqrt(radius**2 - ((r3-r1-dt)**2) / 4)], radius, rpoly1)
            
            da1dr1 = (cura1plusdt-cura1)/dt
            
            r1 = max(r3 - 2*radius, r1 + (a1-cura1)/da1dr1)
            center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
            cura1, _ = getCircleIntersectArea(center, radius, rpoly1)
        
        #adjust radius to match a2
        radius = abs(r3-r1)
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura2, _ = getCircleIntersectArea(center, radius, rpoly2)
        rgap = radius
        doConverge = False
        #higher a = higher cura2
        largeenough = False
        while abs(cura2 - a2)/poly2area > epsilon*(scaleEpsilon**numcycles):
            if not(largeenough):
                largeenough2 = True
                for rpoly2vertex in rpoly2:
                    if rpoly2vertex[1] > 0 and getDistance(rpoly2vertex, center) > radius:
                        largeenough2 = False
                if largeenough2:
                    largeenough = True
            if cura2 > a2:
                radius += rgap
                if not(doConverge):
                    rgap *= 2
                else:
                    rgap /= 2
            else:
                if not(largeenough):
                    radius += rgap
                    rgap *= 2
                elif not(doConverge):
                    rgap /= 4
                    radius -= rgap
                    doConverge = True
                    rgap /= 2
                else:
                    radius -= rgap
                    rgap /= 2
            center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
            cura2, _ = getCircleIntersectArea(center, radius, rpoly2)
            
        #adjust r3 to match a3
        dt = dtbase
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura3, _ = getCircleIntersectArea(center, radius, rpoly3)
        
        while abs(cura3 - a3)/poly3area > epsilon*(scaleEpsilon**numcycles):
            cura3minusdt, _ = getCircleIntersectArea([(r1+r3-dt)/2, math.sqrt(radius**2 - ((r3-dt-r1)**2) / 4)], radius, rpoly3)
            
            da3dr3 = (cura3-cura3minusdt)/dt
            
            r3 = min(r1 + 2*radius, r3 + (a3-cura3)/da3dr3)
            center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
            cura3, _ = getCircleIntersectArea(center, radius, rpoly3)
        
        numcycles += 1
        
        #check for convergence
        cura1, _ = getCircleIntersectArea(center, radius, rpoly1)
        cura2, facetintersects = getCircleIntersectArea(center, radius, rpoly2)
        cura3, _ = getCircleIntersectArea(center, radius, rpoly3)
        if abs(cura1 - a1)/poly1area < epsilon and abs(cura2 - a2)/poly2area < epsilon and abs(cura3 - a3)/poly3area < epsilon:
            returnintersects = [unrot(facetintersect) for _,facetintersect in sorted(zip(list(map(lambda point: point[0], facetintersects)), facetintersects))]
            return unrot(center), radius, returnintersects

    #Max timesteps reached: return failure
    print("Max timesteps reached in getArcFacet({}, {}, {}, {}, {}, {}, {})".format(poly1, poly2, poly3, a1, a2, a3, epsilon))
    return None, None, None

#Newton's method, requires precise initial guess
#Matches area fractions a1, a2, a3
def getArcFacetNewton(poly1, poly2, poly3, a1, a2, a3, epsilon, gcenterx, gcentery, gradius):
    #Center to left of line l1 to l2
    if gradius < 0:
        retcenter, retradius, retintersects = getArcFacetNewton(poly3, poly2, poly1, 1-a3, 1-a2, 1-a1, epsilon, gcenterx, gcentery, -gradius)
        if retcenter is not None:
            retintersects.reverse()
            retradius *= -1
            return retcenter, retradius, retintersects #TODO: radius should always be positive, negative radius means facet is concave, area fraction k -> 1-k
        else:
            return None, None, None

    #Hyperparameters:
    scaleEpsilon = 1 #Adjust threshold by this amount per failed timestep
    dtbase = 1e-8 #Timestep used to calculate numerical estimates of derivatives
    newtonFactor = 1 #Move this amount in Newton direction
    maxTimestep = 500 #Amount of timesteps allowed before we declare failure
    
    #Convert area fractions to areas
    poly1area = getArea(poly1)
    poly2area = getArea(poly2)
    poly3area = getArea(poly3)
    a1 *= poly1area
    a2 *= poly2area
    a3 *= poly3area
    
    #Initial guess
    centerx = gcenterx
    centery = gcentery
    radius = gradius
    
    converged = False
    numcycles = 0

    while (not(converged) and numcycles < maxTimestep):
        #Compute areas
        center = [centerx, centery]
        cura1, _ = getCircleIntersectArea(center, radius, poly1)
        cura2, facetintersects = getCircleIntersectArea(center, radius, poly2)
        cura3, _ = getCircleIntersectArea(center, radius, poly3)
        if abs(cura1 - a1)/poly1area < epsilon and abs(cura2 - a2)/poly2area < epsilon and abs(cura3 - a3)/poly3area < epsilon:
            converged = True
        else:
            dt = dtbase
            
            #Compute derivatives
            xplusdt = [centerx+dt/2, centery]
            cura1plusdt, _ = getCircleIntersectArea(xplusdt, radius, poly1)
            cura2plusdt, _ = getCircleIntersectArea(xplusdt, radius, poly2)
            cura3plusdt, _ = getCircleIntersectArea(xplusdt, radius, poly3)
            xminusdt = [centerx-dt/2, centery]
            cura1minusdt, _ = getCircleIntersectArea(xminusdt, radius, poly1)
            cura2minusdt, _ = getCircleIntersectArea(xminusdt, radius, poly2)
            cura3minusdt, _ = getCircleIntersectArea(xminusdt, radius, poly3)
            
            da1dx = (cura1plusdt-cura1minusdt)/dt
            da2dx = (cura2plusdt-cura2minusdt)/dt
            da3dx = (cura3plusdt-cura3minusdt)/dt
            
            yplusdt = [centerx, centery+dt/2]
            cura1plusdt, _ = getCircleIntersectArea(yplusdt, radius, poly1)
            cura2plusdt, _ = getCircleIntersectArea(yplusdt, radius, poly2)
            cura3plusdt, _ = getCircleIntersectArea(yplusdt, radius, poly3)
            yminusdt = [centerx, centery-dt/2]
            cura1minusdt, _ = getCircleIntersectArea(yminusdt, radius, poly1)
            cura2minusdt, _ = getCircleIntersectArea(yminusdt, radius, poly2)
            cura3minusdt, _ = getCircleIntersectArea(yminusdt, radius, poly3)
            
            da1dy = (cura1plusdt-cura1minusdt)/dt
            da2dy = (cura2plusdt-cura2minusdt)/dt
            da3dy = (cura3plusdt-cura3minusdt)/dt
            
            cura1plusdt, _ = getCircleIntersectArea(center, radius+dt/2, poly1)
            cura2plusdt, _ = getCircleIntersectArea(center, radius+dt/2, poly2)
            cura3plusdt, _ = getCircleIntersectArea(center, radius+dt/2, poly3)
            cura1minusdt, _ = getCircleIntersectArea(center, radius-dt/2, poly1)
            cura2minusdt, _ = getCircleIntersectArea(center, radius-dt/2, poly2)
            cura3minusdt, _ = getCircleIntersectArea(center, radius-dt/2, poly3)
            
            da1dr = (cura1plusdt-cura1minusdt)/dt
            da2dr = (cura2plusdt-cura2minusdt)/dt
            da3dr = (cura3plusdt-cura3minusdt)/dt
            
            jacobian = np.array([[da1dx, da1dy, da1dr], [da2dx, da2dy, da2dr], [da3dx, da3dy, da3dr]])
            det = np.linalg.det(jacobian)
            
            assert det != 0
            
            jacobianinv = np.linalg.inv(jacobian)
            centerx += newtonFactor*(jacobianinv[0][0]*(a1-cura1) + jacobianinv[0][1]*(a2-cura2) + jacobianinv[0][2]*(a3-cura3))
            centery += newtonFactor*(jacobianinv[1][0]*(a1-cura1) + jacobianinv[1][1]*(a2-cura2) + jacobianinv[1][2]*(a3-cura3))
            radius += newtonFactor*(jacobianinv[2][0]*(a1-cura1) + jacobianinv[2][1]*(a2-cura2) + jacobianinv[2][2]*(a3-cura3))

        numcycles += 1
    
    if converged:
        poly1centroid = getCentroid(poly1)
        returnintersects = [facetintersect for _,facetintersect in sorted(zip(list(map(lambda point: getDistance(poly1centroid, point), facetintersects)), facetintersects))]
        return center, radius, returnintersects
    else:
        print("Max timesteps reached in getArcFacetNewton({}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(poly1, poly2, poly3, a1, a2, a3, epsilon, gcenterx, gcentery, gradius))
        return None, None, None