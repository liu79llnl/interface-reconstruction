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
from interface_reconstruction import merge, makeFacets, advectFacets

def advection(gridSize, timesteps, plotVtk, plotMat, makeGapless):
    
    #starting mesh
    #opoints = makeCartesianGrid(gridSize)
    #opoints = makeQuadGrid(gridSize, wiggle)
    #opoints = makeConcaveGrid(gridSize, wiggle)
    opoints = makeFineCartesianGrid(gridSize, resolution)
    
    #mesh velocity per timestep
    velocity = [vsize, vsize]
    
    #variables per time step
    opolys = [[0] * (len(opoints)-1) for _ in range(len(opoints)-1)]

    #matplotlib
    opatches = []
    opatchareas = []
    opatchareapartials = []
    
    #compute initial area setup
    areas = [[0] * len(opolys) for _ in range(len(opolys))]

    for x in range(len(opoints)-1):
        for y in range(len(opoints)-1):
            opoly = [opoints[x][y], opoints[x+1][y], opoints[x+1][y+1], opoints[x][y+1]]
            
            #Initial area setting: "x" and "o" example -> #TODO: to run advection scheme on different shapes, adjust initial area setting
            #Cross
            xpoints = [[4, 0], [10, 6], [16, 0], [20, 4], [14, 10], [20, 16], [16, 20], [10, 14], [4, 20], [0, 16], [6, 10], [0, 4]]
            xpoints = list(map(lambda x : [x[0]+3.005, x[1]+3.025], xpoints))
            xpolyintersects = getPolyIntersectArea(xpoints, opoly)
            for xpolyintersect in xpolyintersects:
                areas[x][y] += abs(getArea(xpolyintersect))
            #Cross
            xpoints = [[6, 0], [14, 0], [14, 6], [20, 6], [20, 14], [14, 14], [14, 20], [6, 20], [6, 14], [0, 14], [0, 6], [6, 6]]
            xpoints = list(map(lambda x : [x[0]+3.005, x[1]+25.025], xpoints))
            xpolyintersects = getPolyIntersectArea(xpoints, opoly)
            for xpolyintersect in xpolyintersects:
                areas[x][y] += abs(getArea(xpolyintersect))
            
            #Union of two circles
            truecircleintersect1, truecircleintersect2 = getCircleCircleIntersects([38.005, 33.005], [38.005, 13.005], 12, 12)
            if getDistance(getCentroid(opoly), truecircleintersect1) <= 2:
                areas[x][y] += getPolyCurvedCornerArea(opoly, [26.005, 33.005], truecircleintersect1, [26.005, 13.005], 12, 12)
            elif getDistance(getCentroid(opoly), truecircleintersect2) <= 2:
                areas[x][y] += getPolyCurvedCornerArea(opoly, [50.005, 13.005], truecircleintersect2, [50.005, 33.005], 12, 12)
            elif getDistance(getCentroid(opoly), [38.005, 33.005]) < getDistance(getCentroid(opoly), [38.005, 13.005]):
                center = [38.005, 33.005]
                area, intersect = getCircleIntersectArea(center, 12, opoly)
                areas[x][y] += area
            else:
                center = [38.005, 13.005]
                area, intersect = getCircleIntersectArea(center, 12, opoly)
                areas[x][y] += area
            
            #Ring
            radiussmall = 7
            center = [38.005, 13.005]
            area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
            areas[x][y] -= area
            #Ring
            radiussmall = 7
            center = [38.005, 33.005]
            area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
            areas[x][y] -= area

            areas[x][y] /= getArea(opoly)
            if areas[x][y] > 1:
                areas[x][y] = 1
            if abs(1-areas[x][y]) < threshold:
                areas[x][y] = 1
            elif abs(areas[x][y]) < threshold:
                areas[x][y] = 0
            elif areas[x][y] < 0:
                areas[x][y] = 0

            opolys[x][y] = opoly

            if plotMat:
                #matplotlib: match partial areas and polygons for plotting
                opatch = Polygon(np.array(opoly), True)
                opatches.append(opatch)
                opatchareas.append(areas[x][y])
                opatchareapartials.append(math.ceil(areas[x][y] - math.floor(areas[x][y])))

    #plot in matplotlib
    if plotMat:
        try:
            os.mkdir('advection_plots')
        except:
            print("Saving plots in ./advection_plots/.")
        op = PatchCollection(opatches, cmap='jet')
        opatchareas = np.array(opatchareas)
        op.set_array(opatchareas)
        fig, ax = plt.subplots()
        ax.set_xlim(-1, len(opolys)/resolution+1)
        ax.set_ylim(-1, len(opolys)/resolution+1)
        ax.add_collection(op)
        plt.savefig("advection_plots/original.png", dpi=199)
        plt.clf()
        
        op = PatchCollection(opatches, cmap='jet')
        opatchareas = np.array(opatchareapartials)
        op.set_array(opatchareas)
        fig, ax = plt.subplots()
        ax.set_xlim(-1, len(opolys)/resolution+1)
        ax.set_ylim(-1, len(opolys)/resolution+1)
        ax.add_collection(op)
        plt.savefig("advection_plots/original_partials.png", dpi=199)
        plt.clf()
    
    #plot in vtk
    if plotVtk:
        try:
            os.mkdir('advection_vtk')
        except:
            print("Saving vtk files in ./advection_vtk/.")
        plotQuadGrid(opoints, 'quads')
    
    #main advection loop
    for timestep in range(timesteps):
        print("Timestep: {}".format(timestep))
        
        mergedpolyindices, mergedpolyinfos, mergedcoords, mergedareafractions = merge(opolys, areas)
        
        predfacets, facetsunique = makeFacets(mergedpolyindices, mergedpolyinfos, mergedcoords, mergedareafractions, opolys, areas, threshold, makeGapless)

        #Plot facets in vtk
        if plotVtk and (timestep % 10 == 0 or timestep == timesteps-1):
            plotFacets(facetsunique, 'timestep_{}'.format(timestep))
        
        print("Computing new areas")
        nareas = advectFacets(opolys, areas, predfacets, velocity, resolution, threshold)
        areas = nareas
        
        if plotMat:
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
            ax.set_xlim(-1, len(opolys)/resolution+1)
            ax.set_ylim(-1, len(opolys)/resolution+1)
            ax.add_collection(newp)
            plt.savefig("advection_plots/pred_timestep_{}.png".format(timestep), dpi=199)
            plt.clf()

            newp = PatchCollection(opatches, cmap='jet')
            npatchareas = []
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    npatchareas.append(math.ceil(areas[x][y] - math.floor(areas[x][y])))
            npatchareas = np.array(npatchareas)
            newp.set_array(npatchareas)
            fig, ax = plt.subplots()
            ax.set_xlim(-1, len(opolys)/resolution+1)
            ax.set_ylim(-1, len(opolys)/resolution+1)
            ax.add_collection(newp)
            plt.savefig("advection_plots/predpartials_timestep_{}.png".format(timestep), dpi=199)
            plt.clf()

            
            #True areas -> #TODO: should correspond to initial area fraction fitting
            trueareas = [[0] * len(opolys) for _ in range(len(opolys))]
            truepatchareas = []
            truepatcherrors = []
            prederrors = []
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    opoly = [opoints[x][y], opoints[x+1][y], opoints[x+1][y+1], opoints[x][y+1]]
                    
                    #Cross
                    xpoints = [[4, 0], [10, 6], [16, 0], [20, 4], [14, 10], [20, 16], [16, 20], [10, 14], [4, 20], [0, 16], [6, 10], [0, 4]]
                    xpoints = list(map(lambda x : [x[0]+3.005+velocity[0]*(timestep+1), x[1]+3.025+velocity[1]*(timestep+1)], xpoints))
                    xpolyintersects = getPolyIntersectArea(xpoints, opoly)
                    for xpolyintersect in xpolyintersects:
                        trueareas[x][y] += abs(getArea(xpolyintersect))
                    #Cross
                    xpoints = [[6, 0], [14, 0], [14, 6], [20, 6], [20, 14], [14, 14], [14, 20], [6, 20], [6, 14], [0, 14], [0, 6], [6, 6]]
                    xpoints = list(map(lambda x : [x[0]+3.005+velocity[0]*(timestep+1), x[1]+25.025+velocity[1]*(timestep+1)], xpoints))
                    xpolyintersects = getPolyIntersectArea(xpoints, opoly)
                    for xpolyintersect in xpolyintersects:
                        trueareas[x][y] += abs(getArea(xpolyintersect))
                    
                    #Union of two circles
                    truecircleintersect1, truecircleintersect2 = getCircleCircleIntersects([38.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)], [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)], 12, 12)
                    if getDistance([x+0.5, y+0.5], truecircleintersect1) <= 2:
                        trueareas[x][y] += getPolyCurvedCornerArea(opoly, [26.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)], truecircleintersect1, [26.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)], 12, 12)
                    elif getDistance([x+0.5, y+0.5], truecircleintersect2) <= 2:
                        trueareas[x][y] += getPolyCurvedCornerArea(opoly, [50.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)], truecircleintersect2, [50.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)], 12, 12)
                    elif getDistance([x+0.5, y+0.5], [38.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)]) < getDistance([x+0.5, y+0.5], [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)]):
                        center = [38.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)]
                        area, intersect = getCircleIntersectArea(center, 12, opoly)
                        trueareas[x][y] += area
                    else:
                        center = [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)]
                        area, intersect = getCircleIntersectArea(center, 12, opoly)
                        trueareas[x][y] += area
                    
                    #Ring
                    radiussmall = 7
                    center = [38.005+velocity[0]*(timestep+1), 13.005+velocity[1]*(timestep+1)]
                    area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
                    trueareas[x][y] -= area
                    #Ring
                    radiussmall = 7
                    center = [38.005+velocity[0]*(timestep+1), 33.005+velocity[1]*(timestep+1)]
                    area, intersect = getCircleIntersectArea(center, radiussmall, opoly)
                    trueareas[x][y] -= area
                    
                    trueareas[x][y] /= getArea(opolys[x][y])
                    if abs(1-trueareas[x][y]) < threshold:
                        trueareas[x][y] = 1
                    elif abs(trueareas[x][y]) < threshold:
                        trueareas[x][y] = 0
                    truepatchareas.append(trueareas[x][y])
                    truepatcherrors.append(math.log10(max(1e-15, abs(trueareas[x][y]-areas[x][y]))))
                    
                    if abs(trueareas[x][y] - areas[x][y]) > 0:
                        prederrors.append(abs(trueareas[x][y] - areas[x][y]))
            
            newp = PatchCollection(opatches, cmap='jet')
            npatchareas = np.array(truepatchareas)
            newp.set_array(npatchareas)
            fig, ax = plt.subplots()
            ax.set_xlim(-1, len(opolys)/resolution+1)
            ax.set_ylim(-1, len(opolys)/resolution+1)
            ax.add_collection(newp)
            plt.savefig("advection_plots/true_timestep_{}.png".format(timestep), dpi=199)
            plt.clf()
            
            plt.hist(list(map(lambda x : math.log10(x), prederrors)))
            plt.xlabel("Predicted area fraction log error")
            plt.ylabel("Frequency")
            plt.title("Errors in predicted area fraction, timestep={}".format(timestep))
            plt.savefig("advection_plots/area_errors_timestep_{}.png".format(timestep))
            plt.clf()

            newp = PatchCollection(opatches, cmap='jet')
            npatchareas = np.array(truepatcherrors)
            newp.set_array(npatchareas)
            fig, ax = plt.subplots()
            ax.set_xlim(-1, len(opolys)/resolution+1)
            ax.set_ylim(-1, len(opolys)/resolution+1)
            ax.add_collection(newp)
            plt.savefig("advection_plots/area_error_patches_{}.png".format(timestep), dpi=199)
            plt.clf()

            prederrors = list(map(lambda x : max(1e-15, x), prederrors))
            prederrors = list(map(lambda x : math.log10(x), prederrors))
            print("Average area fraction error: {}".format(sum(prederrors)/len(prederrors)))

            totaladvectedarea = 0
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    totaladvectedarea += areas[x][y]*getArea(opolys[x][y])

            totaltruearea = 0
            for x in range(len(opolys)):
                for y in range(len(opolys)):
                    totaltruearea += trueareas[x][y]*getArea(opolys[x][y])
                    
            print("Advected area: {}".format(totaladvectedarea))
            print("True area: {}".format(totaltruearea))

            plt.close("all")

np.random.seed(17)
wiggle = 0.5
resolution = 1
threshold = 1e-10*resolution
vsize = 0.0905
gridSize = 70
timesteps = 50
#advection(gridSize, timesteps, plotVtk=True, plotMat=True, makeGapless=True)
advection(gridSize, timesteps, plotVtk=False, plotMat=False, makeGapless=True)