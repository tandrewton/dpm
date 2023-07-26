import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
#from scipy.spatial import Voronoi, voronoi_plot_2d
import pyvoro

# Calculate MSD function
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

# Main function
def main():
    # Set default parameters
    folder = "pipeline/cells/psm/psm_calA01.0_phi0.74_tm10.0_v0"
    v0_arr = ["0.04"]
    k_ecm_arr = ["0.005"]
    k_off_arr = ["1.0"]
    att_arr = ["0.001"]

    for v0 in v0_arr:
        for k_ecm in k_ecm_arr:
            for k_off in k_off_arr:
                for att in att_arr:
                    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
                    folder_path = folder + v0 + "_t_abp1.0k_ecm" + k_ecm + "k_off" + k_off + \
                                  "/_N40_dur1000_att" + att + "_start1_end1_sd1"
                    output_folder = "output" + os.path.sep + folder_path[9:]

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
                        axs[0, 1].plot(xCoords[:, i], yCoords[:, i], '-', linewidth=2)
                    axs[0, 1].set_aspect('equal')
                    axs[0, 1].set_axis_off()
                    #print(np.shape(np.column_stack((xCoords[1, :], yCoords[1, :]))))
                    #print(np.shape([[5,7], [2,3], [3,5],[4,1]]))

                    # Plot Voronoi diagrams
                    for timeii in range(1, 2):  # Change to `xCoords.shape[0]` for all timesteps       
                        cells = pyvoro.compute_2d_voronoi(
                            np.column_stack((xCoords[timeii, :], yCoords[timeii, :])),
                            [[min(xCoords[timeii,:])-1, max(xCoords[timeii,:])+1],[min(yCoords[timeii,:])-1, max(yCoords[timeii,:])+1]], #limits
                            2.0 # block size
                        )
                        print(cells[0])
                    #    vor = Voronoi(np.column_stack((xCoords[timeii, :], yCoords[timeii, :])))
                    #    voronoi_plot_2d(vor, ax=axs[1, 0], show_vertices=False, line_colors='k',
                    #                    line_width=1, line_alpha=0.6)
                    #axs[1, 0].set_aspect('equal')
                    #axs[1, 0].set_axis_off()

                    # want to calculate the neighbor list and compare between frames, then accumulate 

                    # Calculate neighbor distances
                    #neighborDistances = []
                    #for timeii in range(xCoords.shape[0]):
                    #    vor = Voronoi(np.column_stack((xCoords[timeii, :], yCoords[timeii, :])))
                    #    regions, vertices = voronoi_finite_polygons_2d(vor)
                    #    for region in regions:
                    #        indices = region + [region[0]]
                    #        polygon = vertices[indices]
                    #        x, y = np.mean(polygon, axis=0)
                    #        xCoords[timeii].append(x)
                    #        yCoords[timeii].append(y)
                    #        axs[1, 1].plot(polygon[:, 0], polygon[:, 1], 'k-', alpha=0.5)
                    #    axs[1, 1].set_aspect('equal')
                    #    axs[1, 1].set_axis_off()

                    # Save the plot
                    plt.savefig(output_folder + "msd_tracks.png")
                    plt.close()

if __name__ == "__main__":
    main()