import os
import numpy as np
import csv
from io import StringIO
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
from matplotlib.ticker import LinearLocator
from collections import defaultdict
import debugpy

def readCSV(filename):
    # currently omits contacts with the boundary
    with open(filename, 'r') as file:
        csv_content = file.read()
    matrix_strings = csv_content.strip().split('\n\n')
    matrices = []
    for matrix_string in matrix_strings:
        # Create a file-like object to parse each matrix as CSV
        matrix_file = StringIO(matrix_string)

        # Parse the matrix CSV data
        matrix_reader = csv.reader(matrix_file)
        matrix = [list(map(int, row)) for row in matrix_reader]
        matrix = np.array(matrix)
        matrix[matrix > 1] = 1
        matrix = matrix[:-1, :-1] # omit last row, which represents contacts with the boundary

        matrices.append(matrix)
    return matrices

def calculateNeighborExchanges(contactMatrix):
    diff = np.zeros(len(contactMatrix)-1)
    for i in range(len(contactMatrix)-1):
        diff[i] = np.count_nonzero(contactMatrix[i+1] - contactMatrix[i])
        nonzeros = np.nonzero(contactMatrix[i+1] - contactMatrix[i])
        for j in range(len(nonzeros[0])):
            print("frame", i+1, ", nonzero contact indices: ", nonzeros[0][j], nonzeros[1][j])
        debugpy.breakpoint()
    return diff


def main():
    filenames = ["xStream.txt", "xMinStream.txt", "cijStream.txt", "cijMinStream.txt", "comStream.txt"]
    cellPos = np.loadtxt(filenames[0])
    cellMinPos = np.loadtxt(filenames[1])
    #cij = np.loadtxt(filenames[2])
    #cijMin = np.loadtxt(filenames[3])
    cij = readCSV(filenames[2])
    cijMin = readCSV(filenames[3])
    cellCOM = np.loadtxt(filenames[4])

    debugpy.breakpoint()

    # extract header information from cellPos and then remove header
    NCELLS = int(cellPos[0][0])
    NVTOT = int(cellPos[0][1])
    cellPos = cellPos[1:]

    # make movies of cellPos and cellMinPos
    posFig, posAx = plt.subplots()

    # Set up formatting for the movie files
    writervideo = animation.FFMpegWriter(fps=30)

    def updateAnimation(frame, array, rowsPerFrame):
        posAx.clear()
        posAx.set_aspect('equal')  # Set aspect ratio to square
        data = array[frame*rowsPerFrame:(frame+1)*rowsPerFrame]
        posAx.set_xlim(min(data[:, 0]), max(data[:, 0]))  # Set x range
        posAx.set_ylim(min(data[:, 1]), max(data[:, 1]))  # Set y range
        posAx.scatter(data[:,0], data[:,1], s = 80)
        for i in range(NCELLS):
            posAx.text(cellCOM[NCELLS*frame + i][1], cellCOM[NCELLS*frame + i][2], int(cellCOM[NCELLS*frame + i][0]))
        time_text = posAx.text(0.05, 0.95, "")
        time_text.set_text("time = %.1d" % frame)

    debugpy.breakpoint()
    posAni = animation.FuncAnimation(
        posFig, updateAnimation, range(int(cellPos.shape[0]/NVTOT)),fargs=(cellPos, NVTOT), interval=100
    )
    posAni.save('posAnimation.mp4', writer=writervideo, dpi=75)

    minPosAni = animation.FuncAnimation(
        posFig, updateAnimation, range(int(cellPos.shape[0]/NVTOT)), fargs=(cellMinPos, NVTOT), interval=100
    )
    minPosAni.save("minPosAnimation.mp4",  writer=writervideo, dpi=75)

    # use cij and minCij to determine neighbor exchanges as a function of time. minCij should feature fewer switches
    print("calculating NE of positions")
    debugpy.breakpoint()
    NE = calculateNeighborExchanges(cij)
    print("calculating NE of min positions")
    minNE = calculateNeighborExchanges(cijMin)
    figNE, axNE = plt.subplots(1,1)
    axNE.plot(np.arange(0,len(NE)), np.cumsum(NE), 'r')
    axNE.plot(np.arange(0, len(NE)), np.cumsum(minNE), 'k')
    axNE.set_xlabel("Time (min)")
    axNE.set_ylabel("Cumulative NEs")
    plt.show()

    print("closing program")


if __name__ == "__main__":
    main()
