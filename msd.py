import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import pyvoro
from collections import defaultdict
import debugpy


def calculate_msd(filename):
    data = np.loadtxt(filename)
    time = data[:, 0]
    cellPos = data[:, 1:]
    numCells = cellPos.shape[1] // 2
    numSteps = cellPos.shape[0]

    maxLag = numSteps // 4
    msd_per_particle_per_lag = np.zeros((numCells, maxLag))

    for lag in range(1, maxLag+1):
        displacements = cellPos[lag:] - cellPos[:-lag]
        msd_per_particle = np.sum(displacements ** 2, axis=1)
        msd_per_particle_per_lag[:, lag-1] = np.mean(msd_per_particle, axis=0)

    MSD = np.mean(msd_per_particle_per_lag, axis=0)
    return time[1:maxLag+1], MSD


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
                    time, msd = calculate_msd(filename)

                    axs[0, 0].plot(time, msd, 'k')
                    axs[0, 0].set_xlabel("Time (min)")
                    axs[0, 0].set_ylabel("MSD(t)/(25$\pi \mu m^2$)")

                    axs[0, 1].loglog(time, msd, 'k')
                    axs[0, 1].set_xlabel("Time (min)")
                    axs[0, 1].set_ylabel("MSD(t)/(25$\pi \mu m^2$)")
                    axs[0, 1].axline(
                        [10**0, 5*10**-3], [10**1, 5*10**-1], color='r')
                    axs[0, 1].axline(
                        [10**1, 10**-1], [10**2, 10**-0], color='b')
                    # debugpy.breakpoint()

                    # Load data for cell tracks
                    cellPos = np.loadtxt(filename)[:, 1:]
                    xCoords = cellPos[:, ::2]
                    yCoords = cellPos[:, 1::2]

                    for i in range(xCoords.shape[1]):
                        axs[1, 0].plot(
                            xCoords[:, i], yCoords[:, i], '-', linewidth=2)
                    axs[1, 0].set_aspect('equal')
                    axs[1, 0].set_axis_off()

                    # Initialize neighbor, edge data structures for neighbor exchange, edge length, neighbor lifetime analyses
                    neighborExchangeCount = np.zeros(xCoords.shape[0])
                    cellNeighborList = [dict()
                                        for x in range(xCoords.shape[0])]
                    edgeLengthsList = []
                    neighborLifetimesMatrix = np.zeros(
                        (np.shape(xCoords)[1], np.shape(xCoords)[1]))
                    neighborLifetimesList = []
                    deadNeighborSeparationList = []

                    # xCoords.shape[0] = # timesteps
                    for timeii in range(0, xCoords.shape[0]):
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
                                np.arange(0, xCoords.shape[1]-1), list(cellNeighbors[i]))
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
                            # polygon = cell['vertices']
                            # axs[1, 0].fill(
                            #    *zip(*polygon), facecolor='none', edgecolor='k', lw=0.2)
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
                    histLifetime.hist(neighborLifetimesList, bins=2000,
                                      color='blue', alpha=0.7)
                    histLifetime.set_xlabel('Neighbor lifetimes')
                    histLifetime.set_ylabel('Counts')

                    histSeparation = histFig.add_subplot(2, 2, 3)
                    histSeparation.hist(deadNeighborSeparationList, bins=20,
                                        color='blue', alpha=0.7)
                    histSeparation.set_xlabel('Neighbor separations')
                    histSeparation.set_ylabel('Counts')

                    plt.show()

                    persistThreshold = 3  # look 3 configurations forward
                    for timeii in range(0, xCoords.shape[0]):
                        if (timeii > 0 and timeii < xCoords.shape[0] - persistThreshold):
                            # compute differences between neighbors of cell k at time i and time i-1
                            neighborDifferences = {
                                k: cellNeighborList[timeii][k] - cellNeighborList[timeii-1][k] for k in cellNeighborList[timeii]}
                            # for each cell k,
                            for k in neighborDifferences:
                                if (all(element >= 0 for element in neighborDifferences[k])):
                                    # check if neighbor diff nonempty, and k is a new neighbor of some other cell
                                    # note that k is the first cell with a change in neighborlist.
                                    # neighborDifferences[k] is the cell that is a new neighbor of k.
                                    if (neighborDifferences[k] and any(neighborDifferences[l] for l in neighborDifferences[k])):
                                        doesNeighborPersist = True
                                        for l in range(1, persistThreshold):
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
                        range(0, xCoords.shape[0]), np.cumsum(neighborExchangeCount))
                    axs[1, 1].set_xlabel("Time (min)")
                    axs[1, 1].set_ylabel("cumulative NE")

                    # Save the plot
                    plt.show()
                    plt.savefig(output_folder + "msd_tracks.png")
                    plt.close()


if __name__ == "__main__":
    main()
