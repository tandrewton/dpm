import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import debugpy
import csv
from io import StringIO


def readCSV(filename):
    # currently omits contacts with the boundary
    with open(filename, "r") as file:
        csv_content = file.read()
    matrix_strings = csv_content.strip().split("\n\n")
    matrices = []
    for matrix_string in matrix_strings:
        # Create a file-like object to parse each matrix as CSV
        matrix_file = StringIO(matrix_string)

        # Parse the matrix CSV data
        matrix_reader = csv.reader(matrix_file)
        matrix = [list(map(int, row)) for row in matrix_reader]
        matrix = np.array(matrix)
        matrix[matrix > 1] = 1
        # omit last row, which represents contacts with the boundary
        matrix = matrix[:-1, :-1]

        matrices.append(matrix)
    return matrices


def calculateNeighborExchanges(contactMatrix):
    diff = np.zeros(len(contactMatrix) - 1)
    for i in range(len(contactMatrix) - 1):
        diff[i] = np.count_nonzero(contactMatrix[i + 1] - contactMatrix[i])
        nonzeros = np.nonzero(contactMatrix[i + 1] - contactMatrix[i])
    return diff


def filter_df(df, att_arr, att2_arr, v0_arr, tm_arr, gamma_arr):
    # Convert arrays to float and filter the DataFrame
    float_att_arr = [float(i) for i in att_arr]
    float_att2_arr = [float(i) for i in att2_arr]
    float_v0_arr = [float(i) for i in v0_arr]
    float_tm_arr = [float(i) for i in tm_arr]
    float_gamma_arr = [float(i) for i in gamma_arr]

    df_filtered = df[
        df["att"].isin(float_att_arr)
        & df["att2"].isin(float_att2_arr)
        & df["v0"].isin(float_v0_arr)
        & df["tm"].isin(float_tm_arr)
        & df["gamma"].isin(float_gamma_arr)
    ]

    return df_filtered


def process_data(att, v0, att2, tm, gamma, seed, fileheader):
    # for a given simulation (corresponding to a single coordinate in parameter space, and a single seed value)
    #  process the data corresponding to that simulation
    outputFileheader = "output" + fileheader[8:]
    speedFile = outputFileheader + "speed.csv"
    fileExtensions = [
        ".xStream",
        ".xMinStream",
        ".cijStream",
        ".cijMinStream",
        ".comStream",
        ".shapeStream",
    ]
    filenames = [fileheader + ext for ext in fileExtensions]
    # coordinates of all vertices
    try:
        cellPos = np.loadtxt(filenames[0])
        # energy minimized coordinates
        # cellMinPos = np.loadtxt(filenames[1])
        # cij = readCSV(filenames[2])  # cell-cell contact network
        # contact of minimized coordinates
        print(filenames[3])
        cijMin = readCSV(filenames[3])
        # cellCOM = np.loadtxt(filenames[4])
        cellShape = np.loadtxt(filenames[5])
        cellSpeed = np.loadtxt(speedFile)
        assert len(cellSpeed) > 1
    except:
        return pd.DataFrame([]), pd.DataFrame([]), pd.DataFrame([])

    # extract header information from cellPos and then remove header
    # NCELLS = int(cellPos[0][0])
    # NVTOT = int(cellPos[0][1])
    cellPos = cellPos[1:]
    # globalXLim = [min(cellPos[:, 0]), max(cellPos[:, 0])]
    # globalYLim = [min(cellPos[:, 1]), max(cellPos[:, 1])]
    shapeParameters = np.ravel(cellShape[:, 0::2])

    data_shape = {
        "shapeParameter": shapeParameters,
        "att": float(att),
        "v0": float(v0),
        "att2": float(att2),
        "tm": float(tm),
        "gamma": float(gamma),
        "seed": seed,
    }
    df_shape = pd.DataFrame(data_shape)

    minNE = np.mean(calculateNeighborExchanges(cijMin))

    # minNE = 0
    data_NE = {
        "minNE": minNE,
        "att": float(att),
        "v0": float(v0),
        "att2": float(att2),
        "tm": float(tm),
        "gamma": float(gamma),
        "seed": seed,
    }
    df_NE = pd.DataFrame([data_NE])

    # converting speeds to real units since it's not done elsewhere
    data_speed = {
        "speed": cellSpeed * np.sqrt(25 * np.pi) / 3,
        "att": float(att),
        "v0": float(v0),
        "att2": float(att2),
        "tm": float(tm),
        "gamma": float(gamma),
        "seed": seed,
    }
    df_speed = pd.DataFrame(data_speed)

    return df_shape, df_NE, df_speed
    # return df_shape, df_NE
    # return df_shape


def main():
    calA0 = "1.0"
    att_arr = ["0.001", "0.005", "0.01", "0.05", "0.1"]
    kl = "1.0"
    ka = "5.0"
    kb = "0.01"
    v0_arr = ["0.1"]
    att2_arr = ["0.0", "0.001", "0.005", "0.01", "0.05", "0.1"]
    tm_arr = ["10000.0"]
    gamma_arr = ["0"]
    t_abp = "1.0"
    phi = "0.8"
    duration = "500"
    NCELLS = 40
    timeBetweenFrames = 3
    seeds = 25

    from sys import platform
    if platform == "linux" or platform == "linux2":
        # linux
        pass
    elif platform == "darwin":
        folder=f"/Users/AndrewTon/Documents/YalePhD/projects/ZebrafishSegmentation/psm_extracellular_calculation/"
        pipeline_folder=f"pipeline/cells/psm/"
        # OS X
    elif platform == "win32":
        folder=f"C:\\Users\\atata\\projects\\psm_extracellular_calculation\\"
        pipeline_folder=f"pipeline\\cells\\psm\\"

    # psm_calA01.15_phi0.8_tm0_v00.1_t_abp1.0/_N40_dur200_att0.06_att20.012_start1_end10_sd10

    # speed file: C:\Users\atata\projects\dpm\output\cells\psm\psm_calA01.0_phi0.8_tm1.0_v00.0_t_abp1.0_gamma0.5_kl1.0_ka5.0_kb0.1\_N40_dur100_att0.001_att20.001_start1_end1_sd1speed.csv

    # psm_calA01.0_phi0.6_tm10000.0_v00.05_t_abp1.0_gamma0_kl1.0_ka5.0_kb0.1

    df_shapes = pd.DataFrame()
    df_NEs = pd.DataFrame()
    df_speeds = pd.DataFrame()
    # load packing fraction data into dataframe
    df_phis = pd.read_csv(
        f"{folder}windowedPhiDataFrame_calA{calA0}_phi{phi}.txt"
    )

    # load shape and neighbor exchange data into dataframes
    for att in att_arr:
        for v0 in v0_arr:
            for att2 in att2_arr:
                for tm in tm_arr:
                    for gamma in gamma_arr:
                        for seed in range(1, seeds + 1):
                            if platform == "win32":
                                folderStr = f"psm_calA0{calA0}_phi{phi}_tm{tm}_v0{v0}_t_abp{t_abp}_gamma{gamma}_kl{kl}_ka{ka}_kb{kb}\\_N40_dur{duration}_att{att}_att2{att2}_start1_end{seeds}_sd{seed}"
                            else:
                                folderStr = f"psm_calA0{calA0}_phi{phi}_tm{tm}_v0{v0}_t_abp{t_abp}_gamma{gamma}_kl{kl}_ka{ka}_kb{kb}/_N40_dur{duration}_att{att}_att2{att2}_start1_end{seeds}_sd{seed}"
                            fileheader = pipeline_folder + folderStr
                            print(fileheader)

                            # df_shape, df_NE = process_data(att, v0, att2, seed, fileheader)
                            df_shape, df_NE, df_speed = process_data(
                                att, v0, att2, tm, gamma, seed, fileheader
                            )

                            # if dfs returned from process_data are not empty, concat them
                            if (
                                not df_shape.empty
                                and not df_NE.empty
                                and not df_speed.empty
                            ):
                                df_shapes = pd.concat([df_shapes, df_shape])
                                df_NE["minNE"] = (
                                    df_NE["minNE"] / NCELLS / timeBetweenFrames
                                )

                                df_NEs = pd.concat([df_NEs, df_NE])
                                df_speeds = pd.concat([df_speeds, df_speed])
                            else:
                                print(f"{folderStr} is empty!")

    print("finished reading files!")

    # ensure df is restricted to the parameter combinations within the different parameter arrays specified here
    # necessary because the data file may have parameters I'm not interested in
    df_phis = filter_df(df_phis, att_arr, att2_arr, v0_arr, tm_arr, gamma_arr)
    df_phis_grouped_means = (
        df_phis.groupby(["att", "att2", "v0", "tm", "gamma"]).mean().reset_index()
    )

    # Plot all curves on the same axes using hue
    plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
    sns.lineplot(data=df_phis_grouped_means, x="v0", y="phi", marker="o", hue="att")

    print(f"att{att},att2{att2},v0{v0},tm{tm},gamma{gamma},seed{seed}")
    print(df_shapes)

    df_shapes = filter_df(df_shapes, att_arr, att2_arr, v0_arr, tm_arr, gamma_arr)
    df_shapes_grouped_means = (
        df_shapes.groupby(["att", "att2", "v0", "tm", "gamma"]).mean().reset_index()
    )
    # Plot all curves on the same axes using hue
    plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
    sns.lineplot(
        data=df_shapes_grouped_means, x="v0", y="shapeParameter", marker="o", hue="att"
    )

    df_speeds = filter_df(df_speeds, att_arr, att2_arr, v0_arr, tm_arr, gamma_arr)
    df_speeds_grouped_means = (
        df_speeds.groupby(["att", "att2", "v0", "tm", "gamma"]).mean().reset_index()
    )
    # Plot all curves on the same axes using hue
    plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
    sns.lineplot(data=df_speeds_grouped_means, x="v0", y="speed", marker="o", hue="att")

    df_NEs = filter_df(df_NEs, att_arr, att2_arr, v0_arr, tm_arr, gamma_arr)
    df_NEs_grouped_means = (
        df_NEs.groupby(["att", "att2", "v0", "tm", "gamma"]).mean().reset_index()
    )
    # Plot all curves on the same axes using hue
    plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
    sns.lineplot(data=df_NEs_grouped_means, x="v0", y="minNE", marker="o", hue="att")

    # Perform grouping on dataframes in order to make summary plots.

    # Add a column representing unique combinations of parameters.
    # Use lambda here to convert combo column values into strings,
    #  then join the strings into one string so it can be read through seaborn
    # group the dataframe by combo parameters, then average over matching seeds and return a dataframe

    # df_NEs["(att, att2, v0, tm, gamma)"] = df_NEs.apply(lambda row: ', '.join(
    #    [str(row['att']), str(row['att2']), str(row['v0'])]), axis=1)
    # df_NEs_means = df_NEs.groupby(["(att, att2, v0, tm, gamma)", "seed"]).mean().reset_index()

    df_shapes["(att, att2, v0, tm, gamma)"] = df_shapes.apply(
        lambda row: ", ".join(
            [
                str(row["att"]),
                str(row["att2"]),
                str(row["v0"]),
                str(row["tm"]),
                str(row["gamma"]),
            ]
        ),
        axis=1,
    )
    df_shapes_means = (
        df_shapes.groupby(["(att, att2, v0, tm, gamma)", "seed"]).mean().reset_index()
    )

    df_NEs["(att, att2, v0, tm, gamma)"] = df_NEs.apply(
        lambda row: ", ".join(
            [
                str(row["att"]),
                str(row["att2"]),
                str(row["v0"]),
                str(row["tm"]),
                str(row["gamma"]),
            ]
        ),
        axis=1,
    )
    df_NEs_means = (
        df_NEs.groupby(["(att, att2, v0, tm, gamma)", "seed"]).mean().reset_index()
    )

    df_phis["(att, att2, v0, tm, gamma)"] = df_phis.apply(
        lambda row: ", ".join(
            [
                str(row["att"]),
                str(row["att2"]),
                str(row["v0"]),
                str(row["tm"]),
                str(row["gamma"]),
            ]
        ),
        axis=1,
    )
    df_phis_means = (
        df_phis.groupby(["(att, att2, v0, tm, gamma)", "seed"]).mean().reset_index()
    )

    df_speeds["(att, att2, v0, tm, gamma)"] = df_speeds.apply(
        lambda row: ", ".join(
            [
                str(row["att"]),
                str(row["att2"]),
                str(row["v0"]),
                str(row["tm"]),
                str(row["gamma"]),
            ]
        ),
        axis=1,
    )
    df_speeds_means = (
        df_speeds.groupby(["(att, att2, v0, tm, gamma)", "seed"]).mean().reset_index()
    )

    # make boxplots of the mean values for each seed for each parameter combo
    # sns_NEs = sns.catplot(df_NEs_means, x="(att, att2, v0)", y="minNE", kind="box", color="skyblue", height = 8, aspect=2)
    sns_shapes = sns.catplot(
        df_shapes_means,
        x="(att, att2, v0, tm, gamma)",
        y="shapeParameter",
        kind="box",
        color="skyblue",
        height=8,
        aspect=2,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    sns_NEs = sns.catplot(
        df_NEs_means,
        x="(att, att2, v0, tm, gamma)",
        y="minNE",
        kind="box",
        color="skyblue",
        height=8,
        aspect=2,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    sns_phis = sns.catplot(
        df_phis_means,
        x="(att, att2, v0, tm, gamma)",
        y="phi",
        kind="box",
        color="skyblue",
        height=8,
        aspect=2,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    sns_speeds = sns.catplot(
        df_speeds_means,
        x="(att, att2, v0, tm, gamma)",
        y="speed",
        kind="box",
        color="skyblue",
        height=8,
        aspect=2,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    # sns_NEs.figure.savefig(f"sns_NEs_{v0_arr[0]}.png")
    sns_shapes.figure.savefig(f"sns_shapes_{v0_arr[0]}.png")
    sns_NEs.figure.savefig(f"sns_NEs_{v0_arr[0]}.png")
    sns_phis.figure.savefig(f"sns_phis_{v0_arr[0]}.png")
    sns_speeds.figure.savefig(f"sns_speeds_{v0_arr[0]}.png")

    # make boxplots for each parameter combo, pooling observations from each seed
    # sns.catplot(df_NEs, x="(att, att2, v0, tm, gamma)",
    #            y="minNE", kind="box", color="skyblue", height=8, aspect=2)
    sns.catplot(
        df_shapes,
        x="(att, att2, v0, tm, gamma)",
        y="shapeParameter",
        kind="box",
        color="skyblue",
        height=8,
        aspect=2,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    sns.catplot(
        df_NEs,
        x="(att, att2, v0, tm, gamma)",
        y="minNE",
        kind="box",
        color="skyblue",
        height=8,
        aspect=2,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    sns.catplot(
        df_phis,
        x="(att, att2, v0, tm, gamma)",
        y="phi",
        kind="box",
        color="skyblue",
        height=8,
        aspect=2,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    sns.catplot(
        df_speeds,
        x="(att, att2, v0, tm, gamma)",
        y="speed",
        kind="box",
        color="skyblue",
        height=8,
        aspect=2,
    )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    debugpy.breakpoint()

    for fixed_gamma_val in [float(v) for v in gamma_arr]:
        plt.figure()
        fig_heatmap_phis = sns.heatmap(
            df_phis_grouped_means[df_phis_grouped_means["gamma"] == fixed_gamma_val][
                ["att", "att2", "phi"]
            ].pivot(index="att", columns="att2", values="phi"),
            cbar_kws={"label": "$\phi$"},
        )
        plt.xlabel("$\epsilon_2$")
        plt.ylabel("$\epsilon_1$")
        fig_heatmap_phis.invert_yaxis()
        fig_heatmap_phis.figure.savefig(f"fig_heatmap_phis_gamma={fixed_gamma_val}.png")

        plt.figure()
        fig_heatmap_shapes = sns.heatmap(
            df_shapes_grouped_means[
                df_shapes_grouped_means["gamma"] == fixed_gamma_val
            ][["att", "att2", "shapeParameter"]].pivot(
                index="att", columns="att2", values="shapeParameter"
            ),
            cbar_kws={"label": "$\mathcal{A}$"},
        )
        plt.xlabel("$\epsilon_2$")
        plt.ylabel("$\epsilon_1$")
        fig_heatmap_shapes.invert_yaxis()
        fig_heatmap_shapes.figure.savefig(
            f"fig_heatmap_shapes_gamma={fixed_gamma_val}.png"
        )

        plt.figure()
        fig_heatmap_speeds = sns.heatmap(
            df_speeds_grouped_means[
                df_speeds_grouped_means["gamma"] == fixed_gamma_val
            ][["att", "att2", "speed"]].pivot(
                index="att", columns="att2", values="speed"
            ),
            cbar_kws={"label": "v"},
        )
        plt.xlabel("$\epsilon_2$")
        plt.ylabel("$\epsilon_1$")
        fig_heatmap_speeds.invert_yaxis()
        fig_heatmap_speeds.figure.savefig(
            f"fig_heatmap_speeds_gamma={fixed_gamma_val}.png"
        )

        plt.figure()
        fig_heatmap_NEs = sns.heatmap(
            df_NEs_grouped_means[df_NEs_grouped_means["gamma"] == fixed_gamma_val][
                ["att", "att2", "minNE"]
            ].pivot(index="att", columns="att2", values="minNE"),
            cbar_kws={"label": "NE rate"},
        )
        plt.xlabel("$\epsilon_2$")
        plt.ylabel("$\epsilon_1$")
        fig_heatmap_NEs.invert_yaxis()
        fig_heatmap_NEs.figure.savefig(f"fig_heatmap_NEs_gamma={fixed_gamma_val}.png")

    plt.show()


if __name__ == "__main__":
    main()
