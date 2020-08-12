import math
from geoms import getArea, getDistance, getPolyLineArea
from linear_facet import getLinearFacet

#Assumes area is no larger than half the circle
def getArcArea(distance, radius):
    theta = math.asin(distance/(2*radius))
    return (theta * radius**2) - (distance * math.cos(theta)*radius / 2)

#Assumes no duplicate points in a row in poly
def getCircleIntersectArea(center, radius, poly):
    startAt = 1
    apoly = []
    intersects = []
    area = 0
    for i in range(len(poly)):
        p1 = poly[i]
        p2 = poly[(i+1) % len(poly)]
        if getDistance(p1, center) <= radius:
            apoly.append(p1)
        #intersects of circle and segment p1 to p2
        a = (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2
        b = 2*((p1[0]-center[0])*(p2[0]-p1[0]) + (p1[1]-center[1])*(p2[1]-p1[1]))
        c = (p1[0]-center[0])**2 + (p1[1]-center[1])**2 - radius**2
        #a is never 0 because p1 != p2
        disc = b**2 - 4*a*c
        if disc > 0:
            x1 = (-b - math.sqrt(disc))/(2*a)
            if x1 <= 1 and x1 > 0:
                inter1 = [p1[0] + (p2[0]-p1[0])*x1, p1[1] + (p2[1]-p1[1])*x1]
                apoly.append(inter1)
                intersects.append(inter1)
            x2 = (-b + math.sqrt(disc))/(2*a)
            if x2 < 1 and x2 >= 0:
                inter2 = [p1[0] + (p2[0]-p1[0])*x2, p1[1] + (p2[1]-p1[1])*x2]
                apoly.append(inter2)
                intersects.append(inter2)
                if len(intersects) == 1:
                    startAt = 0
    if len(apoly) == 0 and len(intersects) == 0:
        if getDistance(poly[0], center) <= radius:
            return getArea(poly), []
        else:
            return 0, []
    area += getArea(apoly)
    assert len(intersects) % 2 == 0
    for i in range(startAt, len(intersects)+startAt, 2):
        area += getArcArea(getDistance(intersects[i], intersects[(i+1) % len(intersects)]), radius)
    return area, intersects

#Newton's method to match a1, a2, a3s
#Assumes convex polygon
def getArcFacet(poly1, poly2, poly3, a1, a2, a3, epsilon):
    #Hyperparameter:
    scaleEpsilon = 0.99
    dtbase = 1e-8
    
    #Convert area fractions to areas
    poly1area = getArea(poly1)
    poly2area = getArea(poly2)
    poly3area = getArea(poly3)
    a1 *= poly1area
    a2 *= poly2area
    a3 *= poly3area
    
    #Rotate so that x-axis is linear facet
    #print("Getting linear facet:")
    
    l1, l2 = getLinearFacet(poly1, poly3, a1, a3, epsilon/10)
    
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
    
    while not(converged):
        
        #adjust radius to match a2
        radius = abs(r3-r1)
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura2, _ = getCircleIntersectArea(center, radius, rpoly2)
        rgap = radius
        doConverge = False
        #higher a = higher cura2
        while abs(cura2 - a2) > epsilon*(scaleEpsilon**numcycles):
            largeenough = True
            for rpoly2vertex in rpoly2:
                if rpoly2vertex[1] > 0 and getDistance(rpoly2vertex, center) > radius:
                    largeenough = False
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
        
        while abs(cura1 - a1) > epsilon*(scaleEpsilon**numcycles):
            cura1plusdt, _ = getCircleIntersectArea([(r1+dt+r3)/2, math.sqrt(radius**2 - ((r3-r1-dt)**2) / 4)], radius, rpoly1)
            
            da1dr1 = (cura1plusdt-cura1)/dt
            
            r1 = max(r3 - 2*radius, r1 + (a1-cura1)/da1dr1)
            center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
            cura1, _ = getCircleIntersectArea(center, radius, rpoly1)
                
        #adjust r3 to match a3
        dt = dtbase
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura3, _ = getCircleIntersectArea(center, radius, rpoly3)
        
        while abs(cura3 - a3) > epsilon*(scaleEpsilon**numcycles):
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
        if abs(cura1 - a1) < epsilon and abs(cura2 - a2) < epsilon and abs(cura3 - a3) < epsilon:
            returnintersects = []
            for intersect in facetintersects:
                returnintersects.append(unrot(intersect))
            return unrot(center), radius, returnintersects

#No Newton's method
#Assumes convex polygon
def getArcFacet2(poly1, poly2, poly3, a1, a2, a3, epsilon):
    #Hyperparameter:
    scaleEpsilon = 0.99
    
    #Convert area fractions to areas
    poly1area = getArea(poly1)
    poly2area = getArea(poly2)
    poly3area = getArea(poly3)
    a1 *= poly1area
    a2 *= poly2area
    a3 *= poly3area
    
    #Rotate so that x-axis is linear facet
    l1, l2 = getLinearFacet(poly1, poly3, a1, a3, epsilon/10)
    assert getPolyLineArea(poly2, l1, l2) <= a2, "Error in linear facet: min area is {}".format(getPolyLineArea(poly2, l1, l2))

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
    while not(converged):
        
        #adjust radius to match a2
        radius = 1*abs(r3-r1)
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura2, _ = getCircleIntersectArea(center, radius, rpoly2)
        rgap = radius
        doConverge = False
        #higher a = higher cura2
        while abs(cura2 - a2) > epsilon*(scaleEpsilon**numcycles):
            largeenough = True
            for rpoly2vertex in rpoly2:
                if rpoly2vertex[1] > 0 and getDistance(rpoly2vertex, center) > radius:
                    largeenough = False
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
        t1 = 0.5
        r1 = r3 - 2*radius*(1-t1)
        t1gap = 0.25
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura1, _ = getCircleIntersectArea(center, radius, rpoly1)
        doMatcht1 = True
        while doMatcht1:
            if cura1 > a1:
                t1 += t1gap
            else:
                t1 -= t1gap
            t1gap /= 2
            r1 = r3 - 2*radius*(1-t1)
            center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
            cura1, _ = getCircleIntersectArea(center, radius, rpoly1)
            if abs(cura1 - a1) < epsilon/100*(scaleEpsilon**numcycles):
                doMatcht1 = False
                
        #adjust r3 to match a3
        t3 = 0.5
        r3 = r1 + 2*radius*(1-t3)
        t3gap = 0.25
        center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
        cura3, _ = getCircleIntersectArea(center, radius, rpoly3)
        doMatcht3 = True
        while doMatcht3:
            if cura3 > a3:
                t3 += t3gap
            else:
                t3 -= t3gap
            t3gap /= 2
            r3 = r1 + 2*radius*(1-t3)
            center = [(r1+r3)/2, math.sqrt(radius**2 - ((r3-r1)**2) / 4)]
            cura3, _ = getCircleIntersectArea(center, radius, rpoly3)
            if abs(cura3 - a3) < epsilon/100*(scaleEpsilon**numcycles):
                doMatcht3 = False
                
        numcycles += 1
        
        #check for convergence
        cura1, _ = getCircleIntersectArea(center, radius, rpoly1)
        cura2, facetintersects = getCircleIntersectArea(center, radius, rpoly2)
        cura3, _ = getCircleIntersectArea(center, radius, rpoly3)
        if abs(cura1 - a1) < epsilon and abs(cura2 - a2) < epsilon and abs(cura3 - a3) < epsilon:
            returnintersects = []
            for intersect in facetintersects:
                returnintersects.append(unrot(intersect))
            return unrot(center), radius, returnintersects