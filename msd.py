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

                    # Plot cell tracks
                    cellPos = np.loadtxt(filename)[:, 1:]
                    xCoords = cellPos[:, ::2]
                    yCoords = cellPos[:, 1::2]

                    for i in range(xCoords.shape[1]):
                        axs[0, 1].plot(
                            xCoords[:, i], yCoords[:, i], '-', linewidth=2)
                    axs[0, 1].set_aspect('equal')
                    axs[0, 1].set_axis_off()

                    # Plot Voronoi diagrams
                    neighborExchangeCount = np.zeros(xCoords.shape[0])
                    # Change to `xCoords.shape[0]` for all timesteps
                    for timeii in range(0, xCoords.shape[0]):
                        cell_neighbors = defaultdict(set)
                        cells = pyvoro.compute_2d_voronoi(
                            np.column_stack(
                                (xCoords[timeii, :], yCoords[timeii, :])),
                            [[min(xCoords[timeii, :])-1, max(xCoords[timeii, :])+1],
                             [min(yCoords[timeii, :])-1, max(yCoords[timeii, :])+1]],  # limits
                            2.0  # block size
                        )
                        axs[1, 0].set_aspect('equal')
                        axs[1, 0].set_axis_off()
                        for i, cell in enumerate(cells):
                            cell_neighbors[i] = {j['adjacent_cell']
                                                 for j in cells[i]['faces']}
                            polygon = cell['vertices']

                            # axs[1, 0].fill(
                            #    *zip(*polygon), facecolor='none', edgecolor='k', lw=0.2)
                        # calculate number of neighbor exchanges
                        if (timeii > 0):
                            neighborDifferences = {
                                k: cell_neighbors[k] ^ old_cell_neighbors[k] for k in cell_neighbors}
                            if (any(bool(difference) for difference in neighborDifferences.values())):
                                debugpy.breakpoint()
                            for k in neighborDifferences:
                                # if neighbor set of k differs between timesteps
                                # and each difference is symmetric with another neighbor set of neighborDifferences[k]
                                if (neighborDifferences[k] and any(neighborDifferences[l] for l in neighborDifferences[k])):
                                    neighborExchangeCount[timeii] += 1

                        old_cell_neighbors = cell_neighbors
                    # debugpy.breakpoint()

                    # account for double counting of neighbor exchanges by symmetry
                    neighborExchangeCount /= 2
                    plt.plot(
                        range(0, xCoords.shape[0]), np.cumsum(neighborExchangeCount))

                    # Save the plot
                    # plt.show()
                    plt.savefig(output_folder + "msd_tracks.png")
                    plt.close()


if __name__ == "__main__":
    main()
