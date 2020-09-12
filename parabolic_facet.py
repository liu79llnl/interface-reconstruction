import math
import numpy as np
from geoms import getCentroid, getArea, getDistance, lerp, getPolyLineArea, getPolyLineIntersects
from linear_facet import getLinearFacet

#Returns area bounded by parabola and line connecting endpoints of arc
def getParabolaArcArea(a, b, c, x1, x2):
    if x1 <= x2:
        return (x2-x1)/2*(a*(x1**2 + x2**2) + b*(x1 + x2) + 2*c) - (a*(x2**3 - x1**3)/3 + b*(x2**2 - x1**2)/2 + c*(x2 - x1))
    else:
        return -getParabolaArcArea(a, b, c, x2, x1)

#Helper function for getParabolaFacetArea, does not retrieve all parabola-line intersection points
def getParabolaLineIntersects(l1, l2, a, b, c):
    intersects = []
    qa = a*(l2[0]-l1[0])**2
    qb = 2*a*l1[0]*(l2[0]-l1[0]) + b*(l2[0]-l1[0]) + l1[1]-l2[1]
    qc = a*l1[0]**2 + b*l1[0] + c - l1[1]
    if qa == 0 and qb != 0:
        x = -qc/qb
        if x < 1 and x > 0:
            inter = lerp(l1, l2, x)
            intersects.append(inter)
    elif qa != 0:
        disc = qb**2 - 4*qa*qc
        if disc > 0:
            x1 = (-qb - math.sqrt(disc))/(2*qa)
            if x1 <= 1 and x1 > 0:
                inter1 = lerp(l1, l2, x1)
                intersects.append(inter1)
                
            x2 = (-qb + math.sqrt(disc))/(2*qa)
            if x2 < 1 and x2 >= 0:
                inter2 = lerp(l1, l2, x2)
                intersects.append(inter2)
    return intersects

#Returns area above ax^2+bx+c and inside polygon
def getParabolaFacetArea(poly, a, b, c):
    #Hyperparameter:
    adjustcorneramount = 1e-14
    notmod = True
    parabola = lambda x : a*x**2 + b*x + c
    while notmod:
        startAt = 1
        intersectpoints = []
        parabolapoints = []
        for i in range(len(poly)):
            curpoint = poly[i]
            nextpoint = poly[(i+1) % len(poly)]
            curin = curpoint[1] >= parabola(curpoint[0])
            nextin = nextpoint[1] >= parabola(nextpoint[0])
            if curin:
                intersectpoints.append(curpoint)
            lineintersects = getParabolaLineIntersects(curpoint, nextpoint, a, b, c)
            for intersect in lineintersects:
                intersectpoints.append(intersect)
                if len(parabolapoints) == 0 and curin and not(nextin):
                    startAt = 0
                parabolapoints.append(intersect)
        #If not 0 mod 2, parabola intersects a corner, perturb poly and rerun
        if len(parabolapoints) % 2 == 1:
            poly = list(map(lambda x : [x[0]+adjustcorneramount, x[1]+adjustcorneramount], poly))
        else:
            notmod = False

    area = 0
    for i in range(0, len(parabolapoints), 2):
        area += getParabolaArcArea(a, b, c, parabolapoints[startAt+i][0], parabolapoints[(startAt+i+1) % len(parabolapoints)][0])
    area += getArea(intersectpoints)

    return area, parabolapoints

#Matches area fractions a1, a2, a3
#Returns: a, b, c, numpoints parabola facet points, rotation to make parabola x^2 form
#TODO: inner 1D-optimizations can be made a bit quicker via Newton's method
def getParabolaFacet(poly1, poly2, poly3, a1, a2, a3, epsilon, numpoints):
    #Basic sanity tests:
    assert a1 >= 0 and a1 <= 1 and a2 >= 0 and a2 <= 1 and a3 >= 0 and a3 <= 1, print("Given areas for parabola facet are not valid")

    #Hyperparameters:
    scaleEpsilon = 0.99 #Adjust threshold by this amount per failed timestep
    dtbase = 1e-8 #Timestep used to calculate numerical estimates of derivatives
    fixLinearFacetOrientation = epsilon #Threshold used to determine orientation of wrong linear facets
    fixLinearFacet = 1.1 #Amount to adjust linear facet guess by to fix
    maxTimestep = 500 #Amount of timesteps allowed before we declare failure
    
    #Rotate so that x-axis is linear facet
    l1, l2 = getLinearFacet(poly1, poly3, a1, a3, epsilon/10)

    #Convert area fractions to areas
    poly1area = getArea(poly1)
    poly2area = getArea(poly2)
    poly3area = getArea(poly3)
    a1 *= poly1area
    a2 *= poly2area
    a3 *= poly3area

    #Center to left of line l1 to l2
    if getPolyLineArea(poly2, l1, l2) > a2:
        a, b, c, parabolafacet, rot = getParabolaFacet(poly3, poly2, poly1, 1-a3, 1-a2, 1-a1, epsilon, numpoints)
        if a is not None:
            parabolafacet.reverse()
            return -a, b, c, parabolafacet, rot #TODO: a should always be positive, negative a means facet is concave, area fraction k -> 1-k
        else:
            return None, None, None, None, None

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
    
    numcycles = 0
    
    while (not(converged) and numcycles < maxTimestep):
        
        #adjust a to match a2
        amax = float("inf")
        for i in range(len(rpoly2)):
            p1 = rpoly2[i]
            if (p1[0] < r1 or p1[0] > r3) and p1[1] > 0:
                amax = min(amax, p1[1]/((p1[0]-r1)*(p1[0]-r3)))
        #binary search applicable:
        if not(math.isinf(amax)):
            amin = 0
            a = (amin + amax)/2
            agap = (amax - amin)/4
            cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
            #higher a = higher cura2
            while abs(cura2 - a2) > epsilon*(scaleEpsilon**numcycles):
                if cura2 > a2:
                    a -= agap
                else:
                    a += agap
                agap /= 2
                cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
        #binary search not applicable
        else:
            doConverge = False
            a = 0
            agap = 1
            cura2 = getPolyLineArea(poly2, l1, l2)
            while abs(cura2 - a2) > epsilon*(scaleEpsilon**numcycles):
                if cura2 > a2:
                    if not(doConverge):
                        agap /= 4
                        a -= agap
                        doConverge = True
                        agap /= 2
                    else:
                        a -= agap
                        agap /= 2
                else:
                    a += agap
                    if not(doConverge):
                        agap *= 2
                    else:
                        agap /= 2
                cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
            
        #adjust r1 to match a1
        r1down = r3
        r1up = r3
        for i in range(len(rpoly1)):
            p1 = rpoly1[i]
            if p1[0] < r3 and p1[1] < 0:
                r1up = min(r1up, p1[0] - (p1[1]/(a*(p1[0]-r3))))
        t1 = 0.5
        r1 = r1down*t1 + r1up*(1-t1)
        t1gap = 0.25
        cura1, _ = getParabolaFacetArea(rpoly1, a, -a*(r1+r3), a*r1*r3)
        doMatcht1 = True
        while doMatcht1:
            if cura1 > a1:
                t1 += t1gap
            else:
                t1 -= t1gap
            t1gap /= 2
            r1 = r1down*t1 + r1up*(1-t1)
            cura1, _ = getParabolaFacetArea(rpoly1, a, -a*(r1+r3), a*r1*r3)
            if abs(cura1 - a1) < epsilon*(scaleEpsilon**numcycles):
                doMatcht1 = False

        #adjust a to match a2
        amax = float("inf")
        for i in range(len(rpoly2)):
            p1 = rpoly2[i]
            if (p1[0] < r1 or p1[0] > r3) and p1[1] > 0:
                amax = min(amax, p1[1]/((p1[0]-r1)*(p1[0]-r3)))
        #binary search applicable:
        if not(math.isinf(amax)):
            amin = 0
            a = (amin + amax)/2
            agap = (amax - amin)/4
            cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
            #higher a = higher cura2
            while abs(cura2 - a2) > epsilon*(scaleEpsilon**numcycles):
                if cura2 > a2:
                    a -= agap
                else:
                    a += agap
                agap /= 2
                cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
        #binary search not applicable
        else:
            doConverge = False
            a = 0
            agap = 1
            cura2 = getPolyLineArea(poly2, l1, l2)
            while abs(cura2 - a2) > epsilon*(scaleEpsilon**numcycles):
                if cura2 > a2:
                    if not(doConverge):
                        agap /= 4
                        a -= agap
                        doConverge = True
                        agap /= 2
                    else:
                        a -= agap
                        agap /= 2
                else:
                    a += agap
                    if not(doConverge):
                        agap *= 2
                    else:
                        agap /= 2
                cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
            
        #adjust r3 to match a3
        r3down = r1
        r3up = r1
        for i in range(len(rpoly3)):
            p1 = rpoly3[i]
            if p1[0] > r1 and p1[1] < 0:
                r3up = max(r3up, p1[0] - (p1[1]/(a*(p1[0]-r1))))
        t3 = 0.5
        r3 = r3down*t3 + r3up*(1-t3)
        t3gap = 0.25
        cura3, _ = getParabolaFacetArea(rpoly3, a, -a*(r1+r3), a*r1*r3)
        doMatcht3 = True
        while doMatcht3:
            if cura3 > a3:
                t3 += t3gap
            else:
                t3 -= t3gap
            t3gap /= 2
            r3 = r3down*t3 + r3up*(1-t3)
            cura3, _ = getParabolaFacetArea(rpoly3, a, -a*(r1+r3), a*r1*r3)
            if abs(cura3 - a3) < epsilon*(scaleEpsilon**numcycles):
                doMatcht3 = False
        
        numcycles += 1
        
        #check for convergence
        cura1, _ = getParabolaFacetArea(rpoly1, a, -a*(r1+r3), a*r1*r3)
        cura2, facetintersects = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
        cura3, _ = getParabolaFacetArea(rpoly3, a, -a*(r1+r3), a*r1*r3)
        if abs(cura1 - a1) < epsilon and abs(cura2 - a2) < epsilon and abs(cura3 - a3) < epsilon:
            parabolafacet = []
            parabola = lambda x : a*x**2 - a*(r1+r3)*x + a*r1*r3
            for i in range(numpoints):
                parabolafacet.append(unrot(lerp(facetintersects[0], facetintersects[-1], i/numpoints)))
            parabolafacet.append(unrot(facetintersects[-1]))
            return a, -a*(r1+r3), a*r1*r3, parabolafacet, rot

    #Max timesteps reached: return failure
    print("Max timesteps reached in getParabolaFacet({}, {}, {}, {}, {}, {}, {}, {})".format(poly1, poly2, poly3, a1, a2, a3, epsilon, numpoints))
    return None, None, None, None, None

#Newton's method, requires precise initial guess
#TODO: current initial guess is same initial guess as from getParabolaFacet -> can run a few loops of getParabolaFacet to get a better guess
#Matches area fractions a1, a2, a3
def getParabolaFacetNewton(poly1, poly2, poly3, a1, a2, a3, epsilon, numpoints):
    #Basic sanity tests:
    assert a1 >= 0 and a1 <= 1 and a2 >= 0 and a2 <= 1 and a3 >= 0 and a3 <= 1, print("Given areas for parabola facet are not valid")

    #Hyperparameters:
    scaleEpsilon = 0.99 #Adjust threshold by this amount per failed timestep
    dtbase = 1e-8 #Timestep used to calculate numerical estimates of derivatives
    fixLinearFacetOrientation = epsilon #Threshold used to determine orientation of wrong linear facets
    fixLinearFacet = 1.1 #Amount to adjust linear facet guess by to fix
    dtbase = 1e-8 #Timestep used to calculate numerical estimates of derivatives
    newtonFactor = 1 #Move this amount in Newton direction
    maxTimestep = 500 #Amount of timesteps allowed before we declare failure
    
    #Rotate so that x-axis is linear facet
    l1, l2 = getLinearFacet(poly1, poly3, a1, a3, epsilon/10)

    #Convert area fractions to areas
    poly1area = getArea(poly1)
    poly2area = getArea(poly2)
    poly3area = getArea(poly3)
    a1 *= poly1area
    a2 *= poly2area
    a3 *= poly3area

    #Center to left of line l1 to l2
    if getPolyLineArea(poly2, l1, l2) > a2:
        retcenter, retradius, retintersects = getArcFacet(poly3, poly2, poly1, 1-a3, 1-a2, 1-a1, epsilon) #TODO: fix
        if retcenter is not None:
            retintersects.reverse()
            retradius *= -1
            return retcenter, retradius, retintersects
        else:
            return None, None, None

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
    
    numcycles = 0

    #adjust a to match a2
    amax = float("inf")
    for i in range(len(rpoly2)):
        p1 = rpoly2[i]
        if (p1[0] < r1 or p1[0] > r3) and p1[1] > 0:
            amax = min(amax, p1[1]/((p1[0]-r1)*(p1[0]-r3)))
    #binary search applicable:
    if not(math.isinf(amax)):
        amin = 0
        a = (amin + amax)/2
        agap = (amax - amin)/4
        cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
        #higher a = higher cura2
        while abs(cura2 - a2) > epsilon*(scaleEpsilon**numcycles):
            if cura2 > a2:
                a -= agap
            else:
                a += agap
            agap /= 2
            cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)
    #binary search not applicable
    else:
        #print("No binary search")
        doConverge = False
        a = 0
        agap = 1
        cura2 = getPolyLineArea(poly2, l1, l2)
        while abs(cura2 - a2) > epsilon*(scaleEpsilon**numcycles):
            if cura2 > a2:
                if not(doConverge):
                    agap /= 4
                    a -= agap
                    doConverge = True
                    agap /= 2
                else:
                    a -= agap
                    agap /= 2
            else:
                a += agap
                if not(doConverge):
                    agap *= 2
                else:
                    agap /= 2
            cura2, _ = getParabolaFacetArea(rpoly2, a, -a*(r1+r3), a*r1*r3)

    b = -a*(r1+r3)
    c = a*r1*r3

    while (not(converged) and numcycles < maxTimestep):
        #Compute areas
        cura1, _ = getParabolaFacetArea(poly1, a, b, c)
        cura2, facetintersects = getParabolaFacetArea(poly2, a, b, c)
        cura3, _ = getParabolaFacetArea(poly3, a, b, c)
        if abs(cura1 - a1)/poly1area < epsilon and abs(cura2 - a2)/poly2area < epsilon and abs(cura3 - a3)/poly3area < epsilon:
            converged = True
        else:
            dt = dtbase
            
            #Compute derivatives
            aplusdt = a+dt
            cura1plusdt, _ = getParabolaFacetArea(poly1, aplusdt, b, c)
            cura2plusdt, _ = getParabolaFacetArea(poly2, aplusdt, b, c)
            cura3plusdt, _ = getParabolaFacetArea(poly3, aplusdt, b, c)
            aminusdt = a-dt
            cura1minusdt, _ = getParabolaFacetArea(poly1, aminusdt, b, c)
            cura2minusdt, _ = getParabolaFacetArea(poly2, aminusdt, b, c)
            cura3minusdt, _ = getParabolaFacetArea(poly3, aminusdt, b, c)
            
            da1da = (cura1plusdt-cura1minusdt)/dt
            da2da = (cura2plusdt-cura2minusdt)/dt
            da3da = (cura3plusdt-cura3minusdt)/dt
            
            bplusdt = b+dt
            cura1plusdt, _ = getParabolaFacetArea(poly1, a, bplusdt, c)
            cura2plusdt, _ = getParabolaFacetArea(poly2, a, bplusdt, c)
            cura3plusdt, _ = getParabolaFacetArea(poly3, a, bplusdt, c)
            bminusdt = b-dt
            cura1minusdt, _ = getParabolaFacetArea(poly1, a, bminusdt, c)
            cura2minusdt, _ = getParabolaFacetArea(poly2, a, bminusdt, c)
            cura3minusdt, _ = getParabolaFacetArea(poly3, a, bminusdt, c)
            
            da1db = (cura1plusdt-cura1minusdt)/dt
            da2db = (cura2plusdt-cura2minusdt)/dt
            da3db = (cura3plusdt-cura3minusdt)/dt
            
            cplusdt = c+dt
            cura1plusdt, _ = getParabolaFacetArea(poly1, a, b, cplusdt)
            cura2plusdt, _ = getParabolaFacetArea(poly2, a, b, cplusdt)
            cura3plusdt, _ = getParabolaFacetArea(poly3, a, b, cplusdt)
            cminusdt = c-dt
            cura1minusdt, _ = getParabolaFacetArea(poly1, a, b, cminusdt)
            cura2minusdt, _ = getParabolaFacetArea(poly2, a, b, cminusdt)
            cura3minusdt, _ = getParabolaFacetArea(poly3, a, b, cminusdt)
            
            da1dc = (cura1plusdt-cura1minusdt)/dt
            da2dc = (cura2plusdt-cura2minusdt)/dt
            da3dc = (cura3plusdt-cura3minusdt)/dt
            
            jacobian = np.array([[da1da, da1db, da1dc], [da2da, da2db, da2dc], [da3da, da3db, da3dc]])
            det = np.linalg.det(jacobian)
            
            assert det != 0
            
            jacobianinv = np.linalg.inv(jacobian)
            a += newtonFactor*(jacobianinv[0][0]*(a1-cura1) + jacobianinv[0][1]*(a2-cura2) + jacobianinv[0][2]*(a3-cura3))
            b += newtonFactor*(jacobianinv[1][0]*(a1-cura1) + jacobianinv[1][1]*(a2-cura2) + jacobianinv[1][2]*(a3-cura3))
            c += newtonFactor*(jacobianinv[2][0]*(a1-cura1) + jacobianinv[2][1]*(a2-cura2) + jacobianinv[2][2]*(a3-cura3))

        numcycles += 1
    
    if converged:
        poly1centroid = getCentroid(poly1)
        returnintersects = [facetintersect for _,facetintersect in sorted(zip(list(map(lambda point: getDistance(poly1centroid, point), facetintersects)), facetintersects))]
        parabolafacet = []
        parabola = lambda x : a*x**2 - a*(r1+r3)*x + a*r1*r3
        for i in range(numpoints):
            parabolafacet.append(unrot(lerp(facetintersects[0], facetintersects[-1], i/numpoints)))
        parabolafacet.append(unrot(facetintersects[-1]))
        return a, b, c, returnintersects, rot
    else:
        print("Max timesteps reached in getParabolaFacetNewton({}, {}, {}, {}, {}, {}, {}, {})".format(poly1, poly2, poly3, a1, a2, a3, epsilon, numpoints))
        return None, None, None, None, None