import math
from geoms import getArea, getDistance, getPolyLineArea
from linear_facet import getLinearFacet

#Returns area bounded by parabola and line connecting endpoints of arc
def getParabolaArcArea(a, b, c, x1, x2):
    #assert x1 <= x2
    return (x2-x1)/2*(a*(x1**2 + x2**2) + b*(x1 + x2) + 2*c) - (a*(x2**3 - x1**3)/3 + b*(x2**2 - x1**2)/2 + c*(x2 - x1))

#Returns area above ax^2+bx+c and inside polygon
def getParabolaFacetArea(polygon, a, b, c):
    apoly = []
    intersects = []
    area = 0
    parabola = lambda x : a*x**2 + b*x + c
    for i in range(len(polygon)):
        p1 = polygon[i]
        if p1[1] >= parabola(p1[0]):
            apoly.append(p1)
        p2 = polygon[(i+1) % len(polygon)]
        qa = a*(p2[0]-p1[0])**2
        qb = 2*a*p1[0]*(p2[0]-p1[0]) + b*(p2[0]-p1[0]) + p1[1]-p2[1]
        qc = a*p1[0]**2 + b*p1[0] + c - p1[1]
        if qa == 0 and qb != 0:
            x = -qc/qb
            if x < 1 and x > 0:
                inter = [p1[0] + (p2[0]-p1[0])*x, p1[1] + (p2[1]-p1[1])*x]
                apoly.append(inter)
                intersects.append(inter)
        elif qa != 0:
            disc = qb**2 - 4*qa*qc
            if disc > 0:
                x1 = (-qb - math.sqrt(disc))/(2*qa)
                
                if x1 <= 1 and x1 > 0:
                    inter2 = [p1[0] + (p2[0]-p1[0])*x1, p1[1] + (p2[1]-p1[1])*x1]
                    apoly.append(inter2)
                    intersects.append(inter2)
                    
                x2 = (-qb + math.sqrt(disc))/(2*qa)
                if x2 < 1 and x2 >= 0:
                    inter1 = [p1[0] + (p2[0]-p1[0])*x2, p1[1] + (p2[1]-p1[1])*x2]
                    apoly.append(inter1)
                    intersects.append(inter1)
                    
    area += getArea(apoly)
    if a == 0:
        return area
    intersects.sort(key = lambda intersects : intersects[0])
    assert len(intersects) % 2 == 0, intersects
    for i in range(0, len(intersects), 2):
        area += getParabolaArcArea(a, b, c, intersects[i][0], intersects[i+1][0])
    return area, intersects

#Assumes convex polygon
#Returns: a, b, c, numpoints parabola facet points, rotation to make parabola x^2 form
def getParabolaFacet(poly1, poly2, poly3, a1, a2, a3, epsilon, numpoints):
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
    print("Getting linear facet:")
    l1, l2 = getLinearFacet(poly1, poly3, a1, a3, epsilon/10)
    assert getPolyLineArea(poly2, l1, l2) <= a2, "Error in linear facet: min area is {}".format(getPolyLineArea(poly2, l1, l2))
    print("Got linear facet")
    
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
    assert len(intersects1) == 2, "Length of intersects1: {}".format(len(intersects1))
    #towards r1down = cura1 decreases
    r1down = max(intersects1[0], intersects1[1])
    r1up = min(intersects1[0], intersects1[1])
    intersects3 = []
    for i in range(len(rpoly3)):
        p1 = rpoly3[i]
        p2 = rpoly3[(i+1) % len(rpoly3)]
        if ((p1[1] > 0) != (p2[1] > 0)):
            intersects3.append(p1[0] + (p2[0]-p1[0])*p1[1]/(p1[1]-p2[1]))
    assert len(intersects3) == 2, "Length of intersects3: {}".format(len(intersects3))
    #towards r3down = cura3 decreases
    r3down = min(intersects3[0], intersects3[1])
    r3up = max(intersects3[0], intersects3[1])

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
                parabolafacet.append(unrot([facetintersects[0][0]+i/numpoints*(facetintersects[1][0]-facetintersects[0][0]), parabola(facetintersects[0][0]+i/numpoints*(facetintersects[1][0]-facetintersects[0][0]))]))
            parabolafacet.append(unrot(facetintersects[1]))
            return a, -a*(r1+r3), a*r1*r3, parabolafacet, rot