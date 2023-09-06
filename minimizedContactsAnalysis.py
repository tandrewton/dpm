import numpy as np
import csv
from io import StringIO
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection
import debugpy
import argparse

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
        #for j in range(len(nonzeros[0])):
        #    print("frame", i+1, ", nonzero contact indices: ", nonzeros[0][j], nonzeros[1][j])
        debugpy.breakpoint()
    return diff


def main():
    parser = argparse.ArgumentParser(description="Parameter inputs",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-a", "--attraction", help="cell-cell adhesion")
    parser.add_argument("-v0", "--activeVelocity", help="active propulsion coefficient")
    parser.add_argument("-e", "--ecm", help="ecm adhesion")
    args = parser.parse_args()
    att = args.attraction
    v0 = args.activeVelocity
    k_ecm = args.ecm
    # fileheader = "test6"
    fileheader = "pipeline\cells\psm\psm_calA01.0_phi0.74_tm10.0_v0"+v0+"_t_abp1.0k_ecm"+k_ecm+"k_off1.0\_N40_dur100_att"+att+"_start1_end1_sd1"
    outputFileheader = "output"+fileheader[8:]
    fileExtensions = [".xStream", ".xMinStream",
                      ".cijStream", ".cijMinStream", ".comStream"]
    filenames = [fileheader + ext for ext in fileExtensions]
    cellPos = np.loadtxt(filenames[0])
    cellMinPos = np.loadtxt(filenames[1])
    #cij = np.loadtxt(filenames[2])
    #cijMin = np.loadtxt(filenames[3])
    cij = readCSV(filenames[2])
    cijMin = readCSV(filenames[3])
    cellCOM = np.loadtxt(filenames[4])

    # extract header information from cellPos and then remove header
    NCELLS = int(cellPos[0][0])
    NVTOT = int(cellPos[0][1])
    cellPos = cellPos[1:]
    globalXLim = [min(cellPos[:,0]), max(cellPos[:,0])]
    globalYLim = [min(cellPos[:,1]), max(cellPos[:,1])]

    # make movies of cellPos and cellMinPos
    posFig, posAx = plt.subplots()

    # Set up formatting for the movie files
    writervideo = animation.FFMpegWriter(fps=30)

    def updateAnimation(frame, array, rowsPerFrame):
        posAx.clear()
        posAx.set_aspect('equal')  # Set aspect ratio to square
        data = array[frame*rowsPerFrame:(frame+1)*rowsPerFrame]
        posAx.set_xlim(globalXLim[0], globalXLim[1])  # Set x range
        posAx.set_ylim(globalYLim[0], globalYLim[1])  # Set y range
        #posAx.scatter(data[:,0], data[:,1], s = data[:,2])
        circles = [plt.Circle((xi, yi), radius=data[0,2], linewidth=0)
                   for xi, yi in zip(data[:,0], data[:,1])]
        debugpy.breakpoint()
        c = PatchCollection(circles)
        posAx.add_collection(c)
        for i in range(NCELLS):
            posAx.text(cellCOM[NCELLS*frame + i][1], cellCOM[NCELLS*frame + i][2], int(cellCOM[NCELLS*frame + i][0]))
        time_text = posAx.text(0.05, 0.95, "")
        time_text.set_text("time = %.1d" % frame)

    posAni = animation.FuncAnimation(
        posFig, updateAnimation, range(int(cellPos.shape[0]/NVTOT)),fargs=(cellPos, NVTOT), interval=100
    )
    posAni.save(outputFileheader+'posAnimation.mp4', writer=writervideo, dpi=75)

    minPosAni = animation.FuncAnimation(
        posFig, updateAnimation, range(int(cellPos.shape[0]/NVTOT)), fargs=(cellMinPos, NVTOT), interval=100
    )
    minPosAni.save(outputFileheader+"minPosAnimation.mp4",
                   writer=writervideo, dpi=75)

    # use cij and minCij to determine neighbor exchanges as a function of time. minCij should feature fewer switches
    print("calculating NE of positions")
    NE = calculateNeighborExchanges(cij)
    print("calculating NE of min positions")
    minNE = calculateNeighborExchanges(cijMin)
    figNE, axNE = plt.subplots(1,1)
    axNE.plot(np.arange(0,len(NE)), np.cumsum(NE), 'r')
    axNE.plot(np.arange(0, len(NE)), np.cumsum(minNE), 'k')
    axNE.set_xlabel("Time (min)")
    axNE.set_ylabel("Cumulative NEs")
    #plt.show()
    figNE.savefig(outputFileheader+"cumNE.png")
    frameDuration = 0.5
    print("mean NE = ", np.mean(NE)/(frameDuration)/NCELLS)
    print("mean minNE = ", np.mean(minNE)/(frameDuration)/NCELLS)
    with open("NE_overall.txt", "a") as text_file:
        text_file.write(str(att)+","+str(v0)+","+str(k_ecm)+"\n")
        text_file.write("mean NE = "+str(np.mean(NE)/(frameDuration)/NCELLS)+"\n")
        text_file.write("mean minNE = "+str(np.mean(minNE)/(frameDuration)/NCELLS)+"\n\n")

    print("closing program")


if __name__ == "__main__":
    main()
