import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
from matplotlib.ticker import LinearLocator
import pyvoro
from collections import defaultdict
import debugpy


def load_cell_positions(filename):
    data = np.loadtxt(filename)
    time = data[:, 0]
    cell_pos = data[:, 1:]
    return time, cell_pos


def calculate_msd(cell_pos, time, max_lag):
    msd_per_particle_per_lag = np.zeros((cell_pos.shape[1]//2, max_lag))
    for lag in range(1, max_lag+1):
        # displacements = np.diff(cell_pos, axis=0)
        displacements = cell_pos[lag:] - cell_pos[:-lag]
        msd_per_particle = np.sum(displacements ** 2, axis=1)
        msd_per_particle_per_lag[:, lag-1] = np.mean(msd_per_particle, axis=0)
    return time[1:max_lag+1], np.mean(msd_per_particle_per_lag, axis=0)


def create_voronoi_plot(cells, fig, ax):
    # fig, ax = plt.subplots()
    for cell in cells:
        polygon = cell['vertices']
        ax.fill(*zip(*polygon), facecolor='none', edgecolor='k', lw=0.2)
    return None


def main():
    # Set default parameters
    folder = "pipeline/cells/psm/psm_calA01.0_phi0.74_tm10.0_v0"
    v0_arr = ["0.01"]  # , "0.02", "0.04", "0.08", "0.16"]
    k_ecm_arr = ["0.005"]  # , "0.05", "0.5", "5"]
    k_off_arr = ["1.0"]
    att_arr = ["0.001"]  # , "0.01", "0.1"]

    for v0 in v0_arr:
        for k_ecm in k_ecm_arr:
            for k_off in k_off_arr:
                for att in att_arr:
                    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
                    folder_path = folder + v0 + "_t_abp1.0k_ecm" + k_ecm + "k_off" + k_off + \
                        "/_N40_dur1000_att" + att + "_start1_end1_sd1"
                    output_folder = "output" + os.path.sep + folder_path[9:]
                    print(output_folder)

                    filename = folder_path + ".msd"
                    time, cellPos = load_cell_positions(filename)
                    # shorten data so code runs faster
                    time = time[0:len(time)//20]
                    cellPos = cellPos[0:len(cellPos)//20]
                    msdtime, msd = calculate_msd(
                        cellPos, time, len(cellPos)//4)

                    axs[0, 0].plot(msdtime, msd, 'k')
                    axs[0, 0].set_xlabel("Time (min)")
                    axs[0, 0].set_ylabel("MSD(t)/(25$\pi \mu m^2$)")

                    axs[0, 1].loglog(msdtime, msd, 'k')
                    axs[0, 1].set_xlabel("Time (min)")
                    axs[0, 1].set_ylabel("MSD(t)/(25$\pi \mu m^2$)")
                    axs[0, 1].axline(
                        [10**0, 5*10**-3], [10**1, 5*10**-1], color='r')
                    axs[0, 1].axline(
                        [10**1, 10**-1], [10**2, 10**-0], color='b')

                    [xCoords, yCoords] = [cellPos[:, ::2], cellPos[:, 1::2]]

                    for i in range(xCoords.shape[1]):
                        axs[1, 0].plot(
                            xCoords[:, i], yCoords[:, i], '-', linewidth=2)
                    axs[1, 0].set_aspect('equal')
                    axs[1, 0].set_axis_off()

                    # Initialize neighbor, edge data structures for neighbor exchange, edge length, neighbor lifetime analyses
                    neighborExchangeCount = np.zeros(len(cellPos))
                    cellNeighborList = [dict()
                                        for x in range(len(cellPos))]
                    edgeLengthsList = []
                    neighborLifetimesMatrix = np.zeros(
                        (np.shape(xCoords)[1], np.shape(xCoords)[1]))
                    neighborLifetimesList = []
                    deadNeighborSeparationList = []

                    for timeii in range(0, len(cellPos)):
                        cellNeighbors = defaultdict(set)
                        cells = pyvoro.compute_2d_voronoi(
                            np.column_stack(
                                (xCoords[timeii, :], yCoords[timeii, :])),
                            [[min(xCoords[timeii, :])-1, max(xCoords[timeii, :])+1],
                             [min(yCoords[timeii, :])-1, max(yCoords[timeii, :])+1]],  # limits
                            2.0  # block size
                        )
                        # axs[1, 0].set_aspect('equal')
                        # axs[1, 0].set_axis_off()
                        for i, cell in enumerate(cells):
                            # store neighbors of cell i
                            cellNeighbors[i] = {j['adjacent_cell']
                                                for j in cell['faces']}
                            cellNotNeighbors = np.setdiff1d(
                                np.arange(0, xCoords.shape[1]), list(cellNeighbors[i]))
                            # for each neighbor of cell i, record that i-neighbor relationship persists
                            # record if any i-neighbor relationships have ended
                            for adjCellID in cellNeighbors[i]:
                                neighborLifetimesMatrix[i][adjCellID] += 1

                            for nonAdjCellID in cellNotNeighbors:
                                if (nonAdjCellID != i):
                                    # if cell i and nonAdjCellID neighbor lifetime is nonzero, then they used to be neighbors until this timepoint
                                    # record the death of a neighbor relationship: its lifetime, neighbor separation
                                    # reset lifetime to 0
                                    if (neighborLifetimesMatrix[i][nonAdjCellID] > 0):
                                        neighborLifetimesList.append(
                                            neighborLifetimesMatrix[i][nonAdjCellID])
                                        cellPos1 = np.array(
                                            [xCoords[timeii][i], yCoords[timeii][i]])
                                        cellPos2 = np.array(
                                            [xCoords[timeii][nonAdjCellID], yCoords[timeii][nonAdjCellID]])
                                        deadNeighborSeparationList.append(
                                            np.sqrt(np.sum((cellPos1 - cellPos2)**2)))
                                        neighborLifetimesMatrix[i][nonAdjCellID] = 0
                            for j, edge in enumerate(cell['faces']):
                                adjCellID = edge['adjacent_cell']
                                # only choose edges that involve real cells
                                if adjCellID >= 0:
                                    # for cell i and a neighbor, calculate the edge length
                                    edgeVert0 = edge['vertices'][0]
                                    edgeVert1 = edge['vertices'][1]
                                    edgeLength = cell['vertices'][edgeVert0] - \
                                        cell['vertices'][edgeVert1]
                                    edgeLength = np.sqrt(np.sum(edgeLength**2))
                                    edgeLengthsList.append(edgeLength)

                        cellNeighborList[timeii] = cellNeighbors

                    edgeLengthsList = np.asarray(edgeLengthsList)
                    # probably have super long edge lengths to filter, or check if cell ID is negative
                    histFig = plt.figure()
                    histEdgeLength = histFig.add_subplot(2, 2, 1)
                    histEdgeLength.hist(edgeLengthsList, bins=20,
                                        color='blue', alpha=0.7)
                    histEdgeLength.set_xlabel('Edge length')
                    histEdgeLength.set_ylabel('Counts')

                    histLifetime = histFig.add_subplot(2, 2, 2)
                    histLifetime.hist(neighborLifetimesList, bins=100,
                                      color='blue', alpha=0.7, log=True, range=(0, 100))
                    histLifetime.set_xlabel('Neighbor lifetimes')
                    histLifetime.set_ylabel('Counts')

                    histSeparation = histFig.add_subplot(2, 2, 3)
                    histSeparation.hist(deadNeighborSeparationList, bins=20,
                                        color='blue', alpha=0.7)
                    histSeparation.set_xlabel('Neighbor separations')
                    histSeparation.set_ylabel('Counts')
                    # plt.show()
                    # plt.savefig(output_folder + "histograms.png")
                    # plt.close()

                    persistThreshold = 1  # look persistThreshold # configurations forward, minimum is 1
                    neighborCounter = 0
                    neighborCounter2 = 0
                    neighborCounter3 = 0
                    for timeii in range(0, len(cellPos)):
                        if (timeii > 0 and timeii < len(cellPos) - persistThreshold):
                            # compute differences between neighbors of cell k at time i and time i-1
                            neighborDifferences = {
                                k: cellNeighborList[timeii][k] - cellNeighborList[timeii-1][k] for k in cellNeighborList[timeii]}
                            # for each cell k,
                            for k in neighborDifferences:
                                neighborCounter += 1
                                if (all(element >= 0 for element in neighborDifferences[k])):
                                    neighborCounter2 += 1
                                    # check if neighbor diff nonempty, and k is a new neighbor of some other cell
                                    # note that k is the first cell with a change in neighborlist.
                                    # neighborDifferences[k] is the cell that is a new neighbor of k.
                                    # ensures we only consider neighbor exchanges where both cells have new neighbors
                                    if (neighborDifferences[k] and any(neighborDifferences[l] for l in neighborDifferences[k])):
                                        # if (neighborDifferences[k]):
                                        neighborCounter3 += 1
                                        doesNeighborPersist = True
                                        for l in range(1, persistThreshold+1):
                                            # isFutureNeighbor is true if overlap > 0 between sets at different times
                                            isFutureNeighbor = bool(neighborDifferences[
                                                k] & cellNeighborList[timeii + l][k])
                                            doesNeighborPersist = doesNeighborPersist and isFutureNeighbor
                                        if (doesNeighborPersist):
                                            neighborExchangeCount[timeii] += 1
                                            # print(
                                            #    "neighbor exchange between cells", k, ",", neighborDifferences[k])

                    # account for double counting of neighbor exchanges by symmetry
                    neighborExchangeCount /= 2
                    axs[1, 1].plot(
                        range(0, len(cellPos)), np.cumsum(neighborExchangeCount))
                    axs[1, 1].set_xlabel("Time (min)")
                    axs[1, 1].set_ylabel("cumulative NE")

                    print(np.shape(neighborLifetimesList))
                    # Save the plot
                    # plt.show()
                    plt.savefig(output_folder + "msd_tracks.png")
                    # plt.close()

                    tessFig, tessAxs = plt.subplots()

                    def make_animation_step(frame):
                        cells = pyvoro.compute_2d_voronoi(
                            np.column_stack(
                                (xCoords[frame, :], yCoords[frame, :])),
                            [[min(xCoords[frame, :])-1, max(xCoords[frame, :])+1],
                             [min(yCoords[frame, :])-1, max(yCoords[frame, :])+1]],  # limits
                            2.0  # block size
                        )
                        tessAxs.clear()
                        create_voronoi_plot(cells, tessFig, tessAxs)

                    ani2 = animation.FuncAnimation(
                        tessFig, make_animation_step, range(len(time)), interval=100)
                    ani2.save('tessAnimation2.avi', writer='pillow', fps=100)
                    plt.show()


if __name__ == "__main__":
    main()
