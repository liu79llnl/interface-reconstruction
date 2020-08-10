import numpy as np

#Generate concave mesh
def makeCartesianGrid(gridSize):
    print("Making quad mesh")
    points = [[0] * (gridSize+1) for _ in range(gridSize+1)]
    for x in range(gridSize+1):
        for y in range(gridSize+1):
            points[x][y] = [x, y]
    print("Done")
    return points

#Generate quad mesh
def makeQuadGrid(gridSize, wiggle=0.25):
    print("Making quad mesh")
    points = [[0] * (gridSize+1) for _ in range(gridSize+1)]
    for x in range(gridSize+1):
        for y in range(gridSize+1):
            points[x][y] = [x + wiggle*np.random.rand(), y + wiggle*np.random.rand()]
    print("Done")
    return points

#Generate concave mesh
def makeConcaveGrid(gridSize, wiggle):
    print("Making quad mesh")
    points = [[0] * (gridSize+1) for _ in range(gridSize+1)]
    for x in range(gridSize+1):
        for y in range(gridSize+1):
            if (x+y) % 2 == 1:
                points[x][y] = [x-wiggle, y-wiggle]
            else:
                points[x][y] = [x, y]
    print("Done")
    return points