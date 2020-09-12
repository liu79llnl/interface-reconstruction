import numpy as np
import math

#Basic geometric methods

#Implementation of shoelace formula
def getArea(poly):
    if len(poly) < 3:
        return 0
    sum = 0
    for i in range(len(poly)-1):
        sum += poly[i+1][1]*poly[i][0] - poly[i+1][0]*poly[i][1]
    sum += poly[0][1]*poly[-1][0] - poly[0][0]*poly[-1][1]
    return sum/2

#Get Euclidean distance
def getDistance(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

#Get centroid of polygon
def getCentroid(poly):
    assert len(poly) >= 3 and getArea(poly) != 0
    centroid = [0, 0]
    area = 0
    for i in range(len(poly)):
        p1 = poly[i]
        p2 = poly[(i+1) % len(poly)]
        trianglearea = (p2[1]*p1[0] - p2[0]*p1[1])/2
        area += trianglearea
        centroid[0] += (p1[0]+p2[0])/3*trianglearea
        centroid[1] += (p1[1]+p2[1])/3*trianglearea
    centroid[0] /= area
    centroid[1] /= area
    return centroid

#Linear interpolation: if alpha = 0, returns p1; if alpha = 1, returns p2
def lerp(p1, p2, alpha):
    return [p1[0]*(1-alpha)+p2[0]*alpha, p1[1]*(1-alpha)+p2[1]*alpha]

#Returns intersection of segments l1 to l2, p1 to p2, if such an intersection exists
#Returns: intersect, tl, tp
def lineIntersect(l1, l2, p1, p2):
    #two lines intersect once
    try:
        if (p2[0]-p1[0])*(l2[1]-l1[1]) - (p2[1]-p1[1])*(l2[0]-l1[0]) != 0:
            #Not parallel: there is an intersect
            tp = ((p1[1]-l1[1])*(l2[0]-l1[0]) - (p1[0]-l1[0])*(l2[1]-l1[1]))/((p2[0]-p1[0])*(l2[1]-l1[1]) - (p2[1]-p1[1])*(l2[0]-l1[0]))
            tl = ((l1[1]-p1[1])*(p2[0]-p1[0]) - (l1[0]-p1[0])*(p2[1]-p1[1]))/((l2[0]-l1[0])*(p2[1]-p1[1]) - (l2[1]-l1[1])*(p2[0]-p1[0]))
            return lerp(p1, p2, tp), tl, tp
        else:
            #Parallel: call it no intersects
            return None, None, None
            
    except:
        print("Error: {}, {}, {}, {}".format(l1, l2, p1, p2))

#Get intersects of line l1 to l2 with convex polygon
def getPolyLineIntersects(poly, l1, l2):
    assert not (l1[0] == l2[0] and l1[1] == l2[1])
    intersects = []
    distances = []
    if l1[0] == l2[0]:
        for i in range(len(poly)):
            p1 = poly[i]
            p2 = poly[(i+1) % len(poly)]
            if (p1[0] < l1[0] and p2[0] > l1[0]) or (p1[0] > l1[0] and p2[0] < l1[0]):
                t = (l1[0] - p1[0])/(p2[0]-p1[0])
                pinter = [p1[0] + (p2[0]-p1[0])*t, p1[1] + (p2[1]-p1[1])*t]
                intersects.append(pinter)
                distances.append(getDistance(l1, pinter))
    else:
        l = lambda x : l1[1] + (l2[1]-l1[1])*(x-l1[0])/(l2[0]-l1[0])
        for i in range(len(poly)):
            p1 = poly[i]
            p2 = poly[(i+1) % len(poly)]
            if (p1[1] > l(p1[0]) and p2[1] < l(p2[0])) or (p1[1] < l(p1[0]) and p2[1] > l(p2[0])):
                t = (p1[1]-l1[1]-(l2[1]-l1[1])*(p1[0]-l1[0])/(l2[0]-l1[0]))/((l2[1]-l1[1])*(p2[0]-p1[0])/(l2[0]-l1[0]) - (p2[1]-p1[1]))
                pinter = [p1[0] + (p2[0]-p1[0])*t, p1[1] + (p2[1]-p1[1])*t]
                intersects.append(pinter)
                distances.append((l2[0]-l1[0])*(pinter[0]-l1[0]) + (l2[1]-l1[1])*(pinter[1]-l1[1]))
    return [x for _,x in sorted(zip(distances, intersects))]

#Get area to left of line within polygon
def getPolyLineArea(poly, l1, l2):
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
                pinter = [p1[0] + (p2[0]-p1[0])*t, p1[1] + (p2[1]-p1[1])*t]
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
                pinter = [p1[0] + (p2[0]-p1[0])*t, p1[1] + (p2[1]-p1[1])*t]
                intersectRegion.append(pinter)
    return getArea(intersectRegion)

#Points on boundary are not in poly
def pointInPoly(p, poly):
    windNumber = 0
    for i in range(len(poly)):
        vcur = poly[i]
        vnext = poly[(i+1) % len(poly)]
        #upward crossing
        if vcur[1] <= p[1] and vnext[1] > p[1] and (vnext[0]-vcur[0])*(p[1]-vcur[1]) - (vnext[1]-vcur[1])*(p[0]-vcur[0]) > 0:
            windNumber += 1
        #downward crossing
        elif vcur[1] > p[1] and vnext[1] <= p[1] and (vnext[0]-vcur[0])*(p[1]-vcur[1]) - (vnext[1]-vcur[1])*(p[0]-vcur[0]) < 0:
            windNumber -= 1
    
    if windNumber == 0:
        return False
    return True

#TODO: algorithm can be improved from O(n^2) to O(n log n) finding intersection points of pairs of edges
def getPolyIntersectArea(poly1, poly2):
    testedintersects = [[False] * len(poly2) for _ in range(len(poly1))]
    intersections = []
    
    intersectpoints = [[None] * len(poly2) for _ in range(len(poly1))]
    t1s = [[0] * len(poly2) for _ in range(len(poly1))]
    t2s = [[0] * len(poly2) for _ in range(len(poly1))]
    
    #Compute intersects, mark segment pairs that don't intersect
    for index1 in range(len(poly1)):
        for index2 in range(len(poly2)):
            p1cur = poly1[index1]
            p1next = poly1[(index1+1) % len(poly1)]
            p2cur = poly2[index2]
            p2next = poly2[(index2+1) % len(poly2)]
            testintersect, testt1, testt2 = lineIntersect(p1cur, p1next, p2cur, p2next)
            if testintersect is not None and testt1 >= 0 and testt1 <= 1 and testt2 >= 0 and testt2 <= 1:
                intersectpoints[index1][index2] = testintersect
                t1s[index1][index2] = testt1
                t2s[index1][index2] = testt2
            else:
                testedintersects[index1][index2] = True
                t1s[index1][index2] = float("inf")
                t2s[index1][index2] = float("inf")
                
    #Compute intersects
    for index1 in range(len(poly1)):
        for index2 in range(len(poly2)):
            if testedintersects[index1][index2] is False:
                vlist = [intersectpoints[index1][index2]]
                curi1 = index1
                curi2 = index2
                #Compute leftward vector
                v1 = [poly1[(index1+1) % len(poly1)][0]-poly1[index1][0], poly1[(index1+1) % len(poly1)][1]-poly1[index1][1]]
                v2 = [poly2[(index2+1) % len(poly2)][0]-poly2[index2][0], poly2[(index2+1) % len(poly2)][1]-poly2[index2][1]]
                if v1[0]*v2[1]-v1[1]*v2[0] > 0:
                    onOne = False
                    t = t2s[index1][index2]
                else:
                    onOne = True
                    t = t1s[index1][index2]
                
                while vlist[0] != vlist[len(vlist)-1] or len(vlist) < 2:
                    if onOne:
                        mint = float("inf")
                        mini = curi2
                        for i in range(len(poly2)):
                            if t1s[curi1][(curi2+i+1) % len(poly2)] < mint and (t1s[curi1][(curi2+i+1) % len(poly2)] > t or (t1s[curi1][(curi2+i+1) % len(poly2)] == t and i < len(poly2)-1)):
                                mini = (curi2+i+1) % len(poly2)
                                mint = t1s[curi1][mini]
                        if mint == float("inf"):
                            #No swap yet, add next poly1 vertex and continue
                            vlist.append(poly1[(curi1+1) % len(poly1)])
                            t = 0
                            curi1 = (curi1+1) % len(poly1)
                        else:
                            testedintersects[curi1][mini] = True
                            #Swap from poly1 to poly2
                            vlist.append(intersectpoints[curi1][mini])
                            t = t2s[curi1][mini]
                            curi2 = mini
                            onOne = not(onOne)
                    else:
                        mint = float("inf")
                        mini = curi1
                        for i in range(len(poly1)):
                            if t2s[(curi1+i+1) % len(poly1)][curi2] < mint and (t2s[(curi1+i+1) % len(poly1)][curi2] > t or (t2s[(curi1+i+1) % len(poly1)][curi2] == t and i < len(poly1)-1)):
                                mini = (curi1+i+1) % len(poly1)
                                mint = t2s[mini][curi2]
                        if mint == float("inf"):
                            #No swap yet, add next poly2 vertex and continue
                            vlist.append(poly2[(curi2+1) % len(poly2)])
                            t = 0
                            curi2 = (curi2+1) % len(poly2)
                        else:
                            testedintersects[mini][curi2] = True
                            #Swap from poly2 to poly1
                            vlist.append(intersectpoints[mini][curi2])
                            t = t1s[mini][curi2]
                            curi1 = mini
                            onOne = not(onOne)
                               
                intersections.append(vlist)
            else:
                testedintersects[index1][index2] = True

    if len(intersections) > 0:
        return intersections
    else:
        poly1inpoly2 = pointInPoly(poly1[0], poly2)
        poly2inpoly1 = pointInPoly(poly2[0], poly1)
        if poly1inpoly2 and not(poly2inpoly1):
            return [poly1]
        elif poly2inpoly1 and not(poly1inpoly2):
            return [poly2]
        else:
            return []

#Used to merge mesh elements
#Only compares vertices of neighboring mesh polygons
def mergePolys(poly1, poly2):
    #No duplicates
    newpoly = []
    for i in range(len(poly1)):
        vertex1 = poly1[i]
        if vertex1 in poly2:
            index = poly2.index(vertex1)
            if poly1[(i+1) % len(poly1)] == poly2[(index-1) % len(poly2)]:
                for j in range(len(poly2)-2):
                    newpoly.append(poly2[(index+1+j) % len(poly2)])
                for j in range(len(poly1)):
                    newpoly.append(poly1[(i+1+j) % len(poly1)])
                counter = 0
                while counter < len(newpoly):
                    if newpoly[counter] == newpoly[(counter+2) % len(newpoly)]:
                        if counter == len(newpoly)-1:
                            newpoly.pop(0)
                            newpoly.pop(0)
                            counter = max(counter-3, 0)
                        else:
                            newpoly.pop(counter)
                            newpoly.pop(counter)
                            counter = max(counter-1, 0)
                    else:
                        counter += 1
                return newpoly

    print("Failure to merge")
    return None