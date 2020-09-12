import math
from geoms import getArea, getCentroid, getDistance, lerp, getPolyLineArea

#Newton's method
#Assumes centroids are different
#Matches area fractions a1, a2
def getLinearFacet(poly1, poly2, a1, a2, epsilon):
    #Basic sanity tests:
    assert a1 >= 0 and a1 <= 1 and a2 >= 0 and a2 <= 1, "Given areas for linear facet are not valid"

    #Hyperparameters:
    dtbase = 1e-8 #Timestep used to calculate numerical estimates of derivatives
    centroidscale = 1 #Scale initial guesses apart from each other

    #Convert from area fractions to areas
    poly1area = getArea(poly1)
    poly2area = getArea(poly2)

    a1 *= poly1area
    a2 *= poly2area
    
    #Initial guess
    centroid1 = getCentroid(poly1)
    centroid2 = getCentroid(poly2)
    #scale centroids away from polygons so algorithm not affected by sharp corners
    l1 = lerp(centroid1, centroid2, 1-centroidscale)
    l2 = lerp(centroid2, centroid1, 1-centroidscale)
    converged = False
    
    while not(converged):
        #Normal
        normal = [(l1[1]-l2[1])/math.sqrt((l1[1]-l2[1])**2 + (l2[0]-l1[0])**2), (l2[0]-l1[0])/math.sqrt((l1[1]-l2[1])**2 + (l2[0]-l1[0])**2)]
        
        #Compute areas
        cura1 = getPolyLineArea(poly1, l1, l2)
        cura2 = getPolyLineArea(poly2, l1, l2)
        if abs(cura1 - a1)/poly1area < epsilon and abs(cura2 - a2)/poly2area < epsilon:
            converged = True
        else:
        
            #Compute derivatives
            dt = dtbase
            
            l1plusdt = [l1[0]+dt/2*normal[0], l1[1]+dt/2*normal[1]]
            cura1plusdt = getPolyLineArea(poly1, l1plusdt, l2)
            cura2plusdt = getPolyLineArea(poly2, l1plusdt, l2)
            l1minusdt = [l1[0]-dt/2*normal[0], l1[1]-dt/2*normal[1]]
            cura1minusdt = getPolyLineArea(poly1, l1minusdt, l2)
            cura2minusdt = getPolyLineArea(poly2, l1minusdt, l2)

            da1dl1 = (cura1plusdt-cura1minusdt)/dt
            da2dl1 = (cura2plusdt-cura2minusdt)/dt

            l2plusdt = [l2[0]+dt/2*normal[0], l2[1]+dt/2*normal[1]]
            cura1plusdt = getPolyLineArea(poly1, l1, l2plusdt)
            cura2plusdt = getPolyLineArea(poly2, l1, l2plusdt)
            l2minusdt = [l2[0]-dt/2*normal[0], l2[1]-dt/2*normal[1]]
            cura1minusdt = getPolyLineArea(poly1, l1, l2minusdt)
            cura2minusdt = getPolyLineArea(poly2, l1, l2minusdt)

            da1dl2 = (cura1plusdt-cura1minusdt)/dt
            da2dl2 = (cura2plusdt-cura2minusdt)/dt

            det = da1dl1*da2dl2-da2dl1*da1dl2
            
            if det == 0:
                
                prevl1 = l1.copy()
                prevl2 = l2.copy()
                #match a1 by moving l1
                t1 = 0
                t1gap = 1
                l1 = [prevl1[0]+t1*normal[0], prevl1[1]+t1*normal[1]]
                cura1 = getPolyLineArea(poly1, l1, l2)
                cura1dir = (cura1 > a1)
                cura1converged = False
                while abs(cura1 - a1)/poly1area > epsilon:
                    if cura1 < a1:
                        t1 -= t1gap
                    else:
                        t1 += t1gap
                    if cura1converged:
                        t1gap /= 2
                    else:
                        t1gap *= 2
                    l1 = [prevl1[0]+t1*normal[0], prevl1[1]+t1*normal[1]]
                    cura1 = getPolyLineArea(poly1, l1, l2)
                    if not(cura1converged) and (cura1 > a1) != cura1dir:
                        cura1converged = True
                        t1gap /= 4
                
                #match a2 by moving l2
                t2 = 0
                t2gap = 1
                l2 = [prevl2[0]+t2*normal[0], prevl2[1]+t2*normal[1]]
                cura2 = getPolyLineArea(poly2, l1, l2)
                cura2dir = (cura2 > a2)
                cura2converged = False
                while abs(cura2 - a2)/poly2area > epsilon:
                    if cura2 < a2:
                        t2 -= t2gap
                    else:
                        t2 += t2gap
                    if cura2converged:
                        t2gap /= 2
                    else:
                        t2gap *= 2
                    l2 = [prevl2[0]+t2*normal[0], prevl2[1]+t2*normal[1]]
                    cura2 = getPolyLineArea(poly2, l1, l2)
                    if not(cura2converged) and (cura2 > a2) != cura2dir:
                        cura2converged = True
                        t2gap /= 4
            
            else:
                #inverse 2x2 = 1/det * [da2dl2, -da1dl2; -da2dl1, da1dl1]
                t1 = 1/det*(da2dl2*(a1-cura1) - da1dl2*(a2-cura2))
                t2 = 1/det*(-da2dl1*(a1-cura1) + da1dl1*(a2-cura2))

                l1 = [l1[0]+t1*normal[0], l1[1]+t1*normal[1]]
                l2 = [l2[0]+t2*normal[0], l2[1]+t2*normal[1]]
            
    return l1, l2

#Alternate linear facet finder, without Newton's method or numerical derivatives
#Assumes centroids are different
#Matches area fraction a1, a2
def getLinearFacet2(poly1, poly2, a1, a2, epsilon):
    #Basic sanity tests:
    assert a1 >= 0 and a1 <= 1 and a2 >= 0 and a2 <= 1, "Given areas for linear facet are not valid"

    #Hyperparameters:
    scaleEpsilon = 1 #Amount to decrease epsilons by as iterations continue to try to avoid infinite loops
    centroidscale = 5 #Amount to scale centroids away from each other
    cyclethreshold = 1e-4 #2-cycles are detected where distance is greater than cyclethreshold

    #Convert from area fractions to areas
    poly1area = getArea(poly1)
    poly2area = getArea(poly2)

    a1 *= poly1area
    a2 *= poly2area
    
    #Find the pair of vertices in poly1, poly2 such that their distance is minimized
    #centroid1 to centroid2 is closest pair, alsomatch1 to alsomatch2 is second closest or tied
    match1 = poly1[0]
    match2 = poly2[0]
    centroiddistance = getDistance(match1, match2)
    for point1 in poly1:
        for point2 in poly2:
            if getDistance(point1, point2) < centroiddistance:
                match1 = point1
                match2 = point2
                centroiddistance = getDistance(point1, point2)
                
    alsomatch1 = None
    alsomatch2 = None
    alsodistance = None
    for point1 in poly1:
        for point2 in poly2:
            if (alsodistance is None or getDistance(point1, point2) < alsodistance) and point1 != match1 and point2 != match2 and point1 != match2 and point2 != match1:
                alsomatch1 = point1
                alsomatch2 = point2
                alsodistance = getDistance(point1, point2)
    
    #normal points to left of segment centroid1 to centroid2
    if centroiddistance < epsilon:
        #poly1 and poly2 share a vertex within epsilon error
        centroid1 = getCentroid(poly1)
        centroid2 = getCentroid(poly2)
        if alsodistance < epsilon:
            oldnormal = [(centroid1[1]-centroid2[1])/math.sqrt((centroid1[1]-centroid2[1])**2 + (centroid2[0]-centroid1[0])**2), (centroid2[0]-centroid1[0])/math.sqrt((centroid1[1]-centroid2[1])**2 + (centroid2[0]-centroid1[0])**2)]
            normal = [alsomatch2[0]-match1[0], alsomatch2[1]-match1[1]]
            normal2 = [match1[0]-alsomatch1[0], match1[1]-alsomatch1[1]]
            if oldnormal[0]*normal2[0] + oldnormal[1]*normal2[1] > oldnormal[0]*normal[0] + oldnormal[1]*normal[1]:
                normal = normal2
        else:
            normal = [(centroid1[1]-centroid2[1])/math.sqrt((centroid1[1]-centroid2[1])**2 + (centroid2[0]-centroid1[0])**2), (centroid2[0]-centroid1[0])/math.sqrt((centroid1[1]-centroid2[1])**2 + (centroid2[0]-centroid1[0])**2)]
    else:
        #no shared vertex
        centroid1 = getCentroid(poly1)
        centroid2 = getCentroid(poly2)
        normal = [(centroid1[1]-centroid2[1])/math.sqrt((centroid1[1]-centroid2[1])**2 + (centroid2[0]-centroid1[0])**2), (centroid2[0]-centroid1[0])/math.sqrt((centroid1[1]-centroid2[1])**2 + (centroid2[0]-centroid1[0])**2)]
    
    #scale centroids away from polygons so algorithm not affected by sharp corners
    centroid1 = lerp(centroid1, centroid2, 1-centroidscale)
    centroid2 = lerp(centroid2, centroid1, 1-centroidscale)

    t1 = 0
    t2 = 0
    converged = False
    numcycles = 0
    
    #This part used to try to break cycles
    prev2l1 = None
    prev2l2 = None
    prevl1 = None
    prevl2 = None
    
    while not(converged):
        #Idea: first adjust t1 to match a1, then adjust t2 to match a2, then repeat until done
        #higher t1, t3 = lower cura1, cura3 respectively
        l1 = [centroid1[0]+t1*normal[0], centroid1[1]+t1*normal[1]]
        l2 = [centroid2[0]+t2*normal[0], centroid2[1]+t2*normal[1]]
        
        #Cycle
        if prev2l1 is not None and getDistance(prev2l1, l1) < epsilon and getDistance(prevl1, l1) > cyclethreshold and getDistance(prev2l2, l2) < epsilon and getDistance(prevl2, l2) > cyclethreshold:
            #Move the centroids further from each other and from the opposite polygon
            centroidscale *= 2
            centroid1 = lerp(centroid1, centroid2, 1-centroidscale)
            centroid2 = lerp(centroid2, centroid1, 1-centroidscale)
            l1 = centroid1
            l2 = centroid2
        
        prev2l1 = prevl1
        prev2l2 = prevl2
        prevl1 = l1
        prevl2 = l2
        
        cura1 = getPolyLineArea(poly1, l1, l2)
        cura2 = getPolyLineArea(poly2, l1, l2)

        if abs(cura1 - a1)/poly1area < epsilon and abs(cura2 - a2)/poly2area < epsilon:
            converged = True
        else:
            #match a1 by moving l1
            t1 = 0
            t1gap = 1
            l1 = [centroid1[0]+t1*normal[0], centroid1[1]+t1*normal[1]]
            cura1 = getPolyLineArea(poly1, l1, l2)
            cura1dir = (cura1 > a1)
            cura1converged = False
            while abs(cura1 - a1)/poly1area > epsilon*(scaleEpsilon**numcycles):
                if cura1 < a1:
                    t1 -= t1gap
                else:
                    t1 += t1gap
                if cura1converged:
                    t1gap /= 2
                else:
                    t1gap *= 2
                l1 = [centroid1[0]+t1*normal[0], centroid1[1]+t1*normal[1]]
                cura1 = getPolyLineArea(poly1, l1, l2)
                if not(cura1converged) and (cura1 > a1) != cura1dir:
                    cura1converged = True
                    t1gap /= 4

            #match a2 by moving l2
            t2 = 0
            t2gap = 1
            l2 = [centroid2[0]+t2*normal[0], centroid2[1]+t2*normal[1]]
            cura2 = getPolyLineArea(poly2, l1, l2)
            cura2dir = (cura2 > a2)
            cura2converged = False
            while abs(cura2 - a2)/poly2area > epsilon*(scaleEpsilon**numcycles):
                if cura2 < a2:
                    t2 -= t2gap
                else:
                    t2 += t2gap
                if cura2converged:
                    t2gap /= 2
                else:
                    t2gap *= 2
                l2 = [centroid2[0]+t2*normal[0], centroid2[1]+t2*normal[1]]
                cura2 = getPolyLineArea(poly2, l1, l2)
                if not(cura2converged) and (cura2 > a2) != cura2dir:
                    cura2converged = True
                    t2gap /= 4

            numcycles += 1
               
    return l1, l2