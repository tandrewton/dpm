import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
from matplotlib.ticker import LinearLocator
from collections import defaultdict
import debugpy

def main():
    filenames = ["xStream.txt", "xMinStream.txt", "cijStream.txt", "cijMinStream.txt"]
    cellPos = np.loadtxt(filenames[0])
    cellMinPos = np.loadtxt(filenames[1])
    cij = np.loadtxt(filenames[2])
    cijMin = np.loadtxt(filenames[3])
    debugpy.breakpoint()

    # extract header information from cellPos and then remove header
    NCELLS = int(cellPos[0][0])
    NVTOT = int(cellPos[0][1])
    cellPos = cellPos[1:]

    # make movies of cellPos and cellMinPos
    posFig, posAx = plt.subplots()

    def updateAnimation(frame, array, rowsPerFrame):
        posAx.clear()
        posAx.set_aspect('equal')  # Set aspect ratio to square
        posAx.set_xlim(0,7)  # Set x range
        posAx.set_ylim(0,7)  # Set y range
        print(frame, rowsPerFrame, (frame+1)*rowsPerFrame, array.shape)
        data = array[frame*rowsPerFrame:(frame+1)*rowsPerFrame]
        posAx.scatter(data[:,0], data[:,1], s = 80)
        time_text = posAx.text(0.05, 0.95, "")
        time_text.set_text("time = %.1d" % frame)

    posAni = animation.FuncAnimation(
        posFig, updateAnimation, range(int(cellPos.shape[0]/NVTOT)),fargs=(cellPos, NVTOT), interval=100
    )
    posAni.save("posAnimation.gif", writer="pillow", fps=100, dpi=50)

    minPosAni = animation.FuncAnimation(
        posFig, updateAnimation, range(int(cellPos.shape[0]/NVTOT)), fargs=(cellMinPos, NVTOT), interval=100
    )
    minPosAni.save("minPosAnimation.gif", writer="pillow", fps=100, dpi=50)
    debugpy.breakpoint()

    # use cij and minCij to determine neighbor exchanges as a function of time. minCij should feature fewer switches

    print("closing program")


if __name__ == "__main__":
    main()
