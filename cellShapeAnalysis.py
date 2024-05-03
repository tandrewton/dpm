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


def filter_df(
    df, att_arr, att2_arr, v0_arr, kl_arr, gamma_arr, kon_arr, koff_arr, kecm_arr
):
    # Convert arrays to float and filter the DataFrame
    float_att_arr = [float(i) for i in att_arr]
    float_att2_arr = [float(i) for i in att2_arr]
    float_v0_arr = [float(i) for i in v0_arr]
    float_kl_arr = [float(i) for i in kl_arr]
    float_gamma_arr = [float(i) for i in gamma_arr]
    float_kon_arr = [float(i) for i in kon_arr]
    float_koff_arr = [float(i) for i in koff_arr]
    float_kecm_arr = [float(i) for i in kecm_arr]

    df_filtered = df[
        df["att"].isin(float_att_arr)
        & df["att2"].isin(float_att2_arr)
        & df["v0"].isin(float_v0_arr)
        & df["kl"].isin(float_kl_arr)
        & df["gamma"].isin(float_gamma_arr)
        & df["k_on"].isin(float_kon_arr)
        & df["k_off"].isin(float_koff_arr)
        & df["k_ecm"].isin(float_kecm_arr)
    ]

    return df_filtered


def process_data(att, v0, att2, kl, gamma, k_on, k_off, k_ecm, seed, fileheader):
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

    # using 1/shapeParameters = circularity
    data_shape = {
        "circularity": 1 / shapeParameters,
        "att": float(att),
        "v0": float(v0),
        "att2": float(att2),
        "kl": float(kl),
        "gamma": float(gamma),
        "k_on": float(k_on),
        "k_off": float(k_off),
        "k_ecm": float(k_ecm),
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
        "kl": float(kl),
        "gamma": float(gamma),
        "k_on": float(k_on),
        "k_off": float(k_off),
        "k_ecm": float(k_ecm),
        "seed": seed,
    }
    df_NE = pd.DataFrame([data_NE])

    # converting speeds to real units since it's not done elsewhere
    data_speed = {
        "speed": cellSpeed,
        "att": float(att),
        "v0": float(v0),
        "att2": float(att2),
        "kl": float(kl),
        "gamma": float(gamma),
        "k_on": float(k_on),
        "k_off": float(k_off),
        "k_ecm": float(k_ecm),
        "seed": seed,
    }
    df_speed = pd.DataFrame(data_speed)

    return df_shape, df_NE, df_speed
    # return df_shape, df_NE
    # return df_shape


def filter_and_group_data(df, filter_params, group_columns):
    """
    Filter the dataframe based on the specified parameters and then group the data,
    calculating the mean for each group.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        filter_params (list): List of parameters to filter by.
        group_columns (list): Columns to group by.

    Returns:
        pd.DataFrame: Grouped and mean aggregated DataFrame.
    """
    df_filtered = filter_df(df, *filter_params)
    return df_filtered.groupby(group_columns).mean().reset_index()


def create_heatmap(
    data, index_col, columns_col, value_col, x_label, y_label, cbar_label, filename
):
    """
    Create a heatmap from the specified DataFrame.

    Parameters:
        data (DataFrame): The data frame from which to generate the heatmap.
        index_col (str): Column name to use as the rows in the heatmap.
        columns_col (str): Column name to use as the columns in the heatmap.
        value_col (str): Column name to use as the value cells in the heatmap.
        title (str): Title for the heatmap.
        cbar_label (str): Label for the color bar.
        filename (str): File path to save the heatmap image.
    """
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(
        data.pivot(index=index_col, columns=columns_col, values=value_col),
        annot=True,  # Add numbers in each cell
        fmt=".2f",  # Formatting numbers inside the heatmap
        # cmap="viridis",  # Color map
        cbar_kws={"label": cbar_label},
    )
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.gca().invert_yaxis()
    # plt.tight_layout()
    ax.tick_params(axis="y", labelrotation=0)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()  # Close the figure to free memory


def main():
    calA0 = "1.0"
    # att_arr = ["0.001", "0.005", "0.01", "0.03", "0.05", "0.1"]
    att_arr = ["0.01", "0.03"]
    kl_arr = ["0.2"]
    ka = "2.5"
    kb = "0.01"
    v0_arr = ["0.1"]
    # att2_arr = ["0.0005", "0.001", "0.005", "0.01", "0.05", "0.1"]
    att2_arr = ["0", "0.05"]
    tm_arr = ["10000.0"]
    gamma_arr = ["0"]
    kon_arr = ["1.0"]
    koff_arr = ["0.1", "100.0"]
    kecm_arr = att2_arr
    # kecm_arr = ["0.001", "0.005", "0.01", "0.05"]
    # kecm_arr = ["0.05"]
    t_abp = "1.0"
    phi = "0.8"
    duration = "300"
    NCELLS = 40
    # timeBetweenFrames is the time between frames in minutes, which can be calculated from the print intervals in my simulation code. Currently 5 frames
    real_time_unit = 0.5  # minutes, tau_abp in real units
    framesPerTimeUnit = 5
    timeBetweenFrames = framesPerTimeUnit * real_time_unit
    seeds = 10

    font = {"size": 22}

    plt.rc("font", **font)

    from sys import platform

    if platform == "linux" or platform == "linux2":
        # linux
        pass
    elif platform == "darwin":
        folder = f"/Users/AndrewTon/Documents/YalePhD/projects/ZebrafishSegmentation/psm_extracellular_calculation/"
        pipeline_folder = f"pipeline/cells/psm/"
        # OS X
    elif platform == "win32":
        folder = f"C:\\Users\\atata\\projects\\psm_extracellular_calculation\\"
        pipeline_folder = f"pipeline\\cells\\psm\\"

    # psm_calA01.15_phi0.8_tm0_v00.1_t_abp1.0/_N40_dur200_att0.06_att20.012_start1_end10_sd10

    # speed file: C:\Users\atata\projects\dpm\output\cells\psm\psm_calA01.0_phi0.8_tm1.0_v00.0_t_abp1.0_gamma0.5_kl1.0_ka5.0_kb0.1\_N40_dur100_att0.001_att20.001_start1_end1_sd1speed.csv

    # psm_calA01.0_phi0.6_tm10000.0_v00.05_t_abp1.0_gamma0_kl1.0_ka5.0_kb0.1

    for kl in kl_arr:
        # load dfs after kl, because I want to stratify results by kl
        df_shapes = pd.DataFrame()
        df_NEs = pd.DataFrame()
        df_speeds = pd.DataFrame()
        # load packing fraction data into dataframe
        df_phis = pd.read_csv(f"{folder}windowedPhiDataFrame_calA{calA0}_phi{phi}.txt")

        # load shape and neighbor exchange data into dataframes
        for att in att_arr:
            for v0 in v0_arr:
                for att2 in att2_arr:
                    for tm in tm_arr:
                        for gamma in gamma_arr:
                            for k_on in kon_arr:
                                for k_off in koff_arr:
                                    # for k_ecm in kecm_arr:
                                    k_ecm = att2
                                    for seed in range(1, seeds + 1):
                                        if platform == "win32":
                                            folderStr = f"psm_calA0{calA0}_phi{phi}_tm{tm}_v0{v0}_t_abp{t_abp}_gamma{gamma}_k_on_{k_on}_k_off_{k_off}_k_ecm_{k_ecm}_kl{kl}_ka{ka}_kb{kb}\\_N{NCELLS}_dur{duration}_att{att}_att2{att2}_start1_end{seeds}_sd{seed}"
                                        else:
                                            folderStr = f"psm_calA0{calA0}_phi{phi}_tm{tm}_v0{v0}_t_abp{t_abp}_gamma{gamma}_k_on_{k_on}_k_off_{k_off}_k_ecm_{k_ecm}_kl{kl}_ka{ka}_kb{kb}/_N{NCELLS}_dur{duration}_att{att}_att2{att2}_start1_end{seeds}_sd{seed}"
                                        fileheader = pipeline_folder + folderStr
                                        print(fileheader)

                                        # df_shape, df_NE = process_data(att, v0, att2, seed, fileheader)
                                        df_shape, df_NE, df_speed = process_data(
                                            att,
                                            v0,
                                            att2,
                                            kl,
                                            gamma,
                                            k_on,
                                            k_off,
                                            k_ecm,
                                            seed,
                                            fileheader,
                                        )

                                        # if dfs returned from process_data are not empty, concat them
                                        if (
                                            not df_shape.empty
                                            and not df_NE.empty
                                            and not df_speed.empty
                                        ):
                                            df_shapes = pd.concat([df_shapes, df_shape])
                                            df_NE["minNE"] = (
                                                df_NE["minNE"]
                                                / NCELLS
                                                / timeBetweenFrames
                                            )

                                            df_speed["speed"] = (
                                                df_speed["speed"]
                                                * np.sqrt(25 * np.pi)
                                                / real_time_unit
                                            )

                                            df_NEs = pd.concat([df_NEs, df_NE])
                                            df_speeds = pd.concat([df_speeds, df_speed])
                                        else:
                                            print(f"{folderStr} is empty!")

        print("finished reading files!")

        # ensure df is restricted to the parameter combinations within the different parameter arrays specified here
        # necessary because the data file may have parameters I'm not interested in
        grouping_list = ["att", "att2", "v0", "kl", "gamma", "k_on", "k_off", "k_ecm"]
        df_phis = filter_df(
            df_phis,
            att_arr,
            att2_arr,
            v0_arr,
            [kl],
            gamma_arr,
            kon_arr,
            koff_arr,
            kecm_arr,
        )
        df_phis_grouped_means = df_phis.groupby(grouping_list).mean().reset_index()
        plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
        sns.lineplot(data=df_phis_grouped_means, x="v0", y="phi", marker="o", hue="att")

        df_shapes = filter_df(
            df_shapes,
            att_arr,
            att2_arr,
            v0_arr,
            [kl],
            gamma_arr,
            kon_arr,
            koff_arr,
            kecm_arr,
        )
        df_shapes_grouped_means = df_shapes.groupby(grouping_list).mean().reset_index()
        plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
        sns.lineplot(
            data=df_shapes_grouped_means, x="v0", y="circularity", marker="o", hue="att"
        )

        df_speeds = filter_df(
            df_speeds,
            att_arr,
            att2_arr,
            v0_arr,
            [kl],
            gamma_arr,
            kon_arr,
            koff_arr,
            kecm_arr,
        )
        df_speeds_grouped_means = df_speeds.groupby(grouping_list).mean().reset_index()
        plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
        sns.lineplot(
            data=df_speeds_grouped_means, x="v0", y="speed", marker="o", hue="att"
        )

        df_NEs = filter_df(
            df_NEs,
            att_arr,
            att2_arr,
            v0_arr,
            [kl],
            gamma_arr,
            kon_arr,
            koff_arr,
            kecm_arr,
        )
        df_NEs_grouped_means = df_NEs.groupby(grouping_list).mean().reset_index()
        plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
        sns.lineplot(
            data=df_NEs_grouped_means, x="v0", y="minNE", marker="o", hue="att"
        )

        # Perform grouping on dataframes in order to make summary plots.

        # Add a column representing unique combinations of parameters.
        # Use lambda here to convert combo column values into strings,
        #  then join the strings into one string so it can be read through seaborn
        # group the dataframe by combo parameters, then average over matching seeds and return a dataframe

        grouping_str = "(att, att2, v0, kl, gamma, k_on, k_off, k_ecm)"
        df_shapes[grouping_str] = df_shapes.apply(
            lambda row: ", ".join(
                [
                    str(row["att"]),
                    str(row["att2"]),
                    str(row["v0"]),
                    str(row["kl"]),
                    str(row["gamma"]),
                    str(row["k_on"]),
                    str(row["k_off"]),
                    str(row["k_ecm"]),
                ]
            ),
            axis=1,
        )
        df_shapes_means = df_shapes.groupby([grouping_str, "seed"]).mean().reset_index()

        df_NEs[grouping_str] = df_NEs.apply(
            lambda row: ", ".join(
                [
                    str(row["att"]),
                    str(row["att2"]),
                    str(row["v0"]),
                    str(row["kl"]),
                    str(row["gamma"]),
                    str(row["k_on"]),
                    str(row["k_off"]),
                    str(row["k_ecm"]),
                ]
            ),
            axis=1,
        )
        df_NEs_means = df_NEs.groupby([grouping_str, "seed"]).mean().reset_index()

        df_phis[grouping_str] = df_phis.apply(
            lambda row: ", ".join(
                [
                    str(row["att"]),
                    str(row["att2"]),
                    str(row["v0"]),
                    str(row["kl"]),
                    str(row["gamma"]),
                    str(row["k_on"]),
                    str(row["k_off"]),
                    str(row["k_ecm"]),
                ]
            ),
            axis=1,
        )
        df_phis_means = df_phis.groupby([grouping_str, "seed"]).mean().reset_index()

        df_speeds[grouping_str] = df_speeds.apply(
            lambda row: ", ".join(
                [
                    str(row["att"]),
                    str(row["att2"]),
                    str(row["v0"]),
                    str(row["kl"]),
                    str(row["gamma"]),
                    str(row["k_on"]),
                    str(row["k_off"]),
                    str(row["k_ecm"]),
                ]
            ),
            axis=1,
        )
        df_speeds_means = df_speeds.groupby([grouping_str, "seed"]).mean().reset_index()

        # make boxplots of the mean values for each seed for each parameter combo
        # sns_NEs = sns.catplot(df_NEs_means, x="(att, att2, v0)", y="minNE", kind="box", color="skyblue", height = 8, aspect=2)
        sns_shapes = sns.catplot(
            df_shapes_means,
            x=grouping_str,
            y="circularity",
            kind="box",
            color="skyblue",
            height=8,
            aspect=2,
        )

        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        sns_NEs = sns.catplot(
            df_NEs_means,
            x=grouping_str,
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
            x=grouping_str,
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
            x=grouping_str,
            y="speed",
            kind="box",
            color="skyblue",
            height=8,
            aspect=2,
        )

        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        # sns_NEs.figure.savefig(f"sns_NEs_{v0_arr[0]}.png")
        sns_shapes.figure.savefig(f"sns_shapes_{v0_arr[0]}_kl{kl}.png")
        sns_NEs.figure.savefig(f"sns_NEs_{v0_arr[0]}_kl{kl}.png")
        sns_phis.figure.savefig(f"sns_phis_{v0_arr[0]}_kl{kl}.png")
        sns_speeds.figure.savefig(f"sns_speeds_{v0_arr[0]}_kl{kl}.png")

        # make boxplots for each parameter combo, pooling observations from each seed
        # sns.catplot(df_NEs, x="(att, att2, v0, kl, gamma)",
        #            y="minNE", kind="box", color="skyblue", height=8, aspect=2)
        sns.catplot(
            df_shapes,
            x=grouping_str,
            y="circularity",
            kind="box",
            color="skyblue",
            height=8,
            aspect=2,
        )

        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        sns.catplot(
            df_NEs,
            x=grouping_str,
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
            x=grouping_str,
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
            x=grouping_str,
            y="speed",
            kind="box",
            color="skyblue",
            height=8,
            aspect=2,
        )

        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        plt.rcParams.update({"font.size": 30})
        # for fixed_val in [float(v) for v in gamma_arr]:
        for fixed_val in [float(v) for v in koff_arr]:
            print("koff = ", fixed_val)
            create_heatmap(
                df_phis_grouped_means[
                    (df_phis_grouped_means["k_off"] == fixed_val)
                    & (df_phis_grouped_means["kl"] == float(kl))
                ][["att", "att2", "v0", "k_ecm", "phi"]],
                "att",
                "att2",
                "phi",
                "$\epsilon_1$",
                "$\epsilon_2$",
                "$\phi$",
                f"fig_heatmap_phis_k_off={fixed_val}_kl={kl}.png",
            )

            create_heatmap(
                df_shapes_grouped_means[
                    (df_shapes_grouped_means["k_off"] == fixed_val)
                    & (df_shapes_grouped_means["kl"] == float(kl))
                ][["att", "att2", "k_ecm", "circularity"]],
                "att",
                "att2",
                "circularity",
                "$\epsilon_1$",
                "$\epsilon_2$",
                "$C$",
                f"fig_heatmap_shapes_k_off={fixed_val}_kl={kl}.png",
            )

            create_heatmap(
                df_speeds_grouped_means[
                    (df_speeds_grouped_means["k_off"] == fixed_val)
                    & (df_speeds_grouped_means["kl"] == float(kl))
                ][["att", "att2", "k_ecm", "speed"]],
                "att",
                "att2",
                "speed",
                "$\epsilon_1$",
                "$\epsilon_2$",
                "$v$",
                f"fig_heatmap_speeds_k_off={fixed_val}_kl={kl}.png",
            )

            create_heatmap(
                df_NEs_grouped_means[
                    (df_NEs_grouped_means["k_off"] == fixed_val)
                    & (df_NEs_grouped_means["kl"] == float(kl))
                ][["att", "att2", "k_ecm", "minNE"]],
                "att",
                "att2",
                "minNE",
                "$\epsilon_1$",
                "$\epsilon_2$",
                "NE rate",
                f"fig_heatmap_NEs_k_off={fixed_val}_kl={kl}.png",
            )

            # make sure to change these parameters, and any of the bracketed f-string components
            #  if those parameters change in the simulations, or if I want to explore different phenotype values
            genotype_tags = []
            for att_val in att_arr[::-1]:
                for att2_val in att2_arr[::-1]:
                    genotype_tags.append(
                        f"{att_val}, {float(att2_val)}, {v0}, {kl}, 0.0, 1.0, {fixed_val}, {k_ecm}"
                    )

            df_all = pd.DataFrame(
                columns=[
                    "params",
                    "phi_mean",
                    "phi_std",
                    "shape_mean",
                    "shape_std",
                    "v_mean",
                    "v_std",
                    "NE_mean",
                    "NE_std",
                ]
            )

            for i in range(0, len(genotype_tags)):
                phi_group = df_phis[
                    df_phis["(att, att2, v0, kl, gamma, k_on, k_off, k_ecm)"]
                    == genotype_tags[i]
                ].groupby("seed")["phi"]
                shape_group = df_shapes[
                    df_shapes["(att, att2, v0, kl, gamma, k_on, k_off, k_ecm)"]
                    == genotype_tags[i]
                ].groupby("seed")["circularity"]
                speed_group = df_speeds[
                    df_speeds["(att, att2, v0, kl, gamma, k_on, k_off, k_ecm)"]
                    == genotype_tags[i]
                ].groupby("seed")["speed"]
                NE_group = df_NEs[
                    df_NEs["(att, att2, v0, kl, gamma, k_on, k_off, k_ecm)"]
                    == genotype_tags[i]
                ].groupby("seed")["minNE"]
                new_row = {
                    "params": genotype_tags[i],
                    "phi_mean": phi_group.mean(),
                    "phi_std": phi_group.std(),
                    "shape_mean": shape_group.mean(),
                    "shape_std": shape_group.std(),
                    "v_mean": speed_group.mean(),
                    "v_std": speed_group.std(),
                    "NE_mean": NE_group.mean(),
                    "NE_std": NE_group.std(),
                }
                df_all = pd.concat([df_all, pd.DataFrame(new_row)], ignore_index=True)

            # Get the unique values in the order they appear in the DataFrame
            ordered_categories = df_all["params"].drop_duplicates().tolist()

            # Convert the 'params' column to a categorical type with the specified order
            df_all["params"] = pd.Categorical(
                df_all["params"], categories=ordered_categories, ordered=True
            )

            # Convert 'Category' to a categorical type and assign numerical values
            df_all["eps1,eps2"] = df_all["params"].astype("category").cat.codes
            offset_width = 0.1  # Adjust the spacing between points in the same category
            grouped = df_all.groupby("eps1,eps2")
            offsets = {
                cat: np.linspace(-offset_width / 2, offset_width / 2, num=len(group))
                for cat, group in grouped
            }
            df_all["Offset"] = df_all.apply(
                lambda row: offsets[row["eps1,eps2"]][
                    grouped.groups[row["eps1,eps2"]].get_loc(row.name)
                ],
                axis=1,
            )

            if fixed_val == float(koff_arr[0]):
                fig, axs = plt.subplots(4, 1, sharex=True, figsize=(15, 15))
                clr = "red"
                constant_offset = 0.2
                # constant_offset = 0
            elif fixed_val == float(koff_arr[1]):
                # comment out this guy below to put figures on same plot for koff_arr[i] i > 0
                # fig, axs = plt.subplots(4, 1, sharex=True, figsize=(15, 15))
                clr = "black"
                constant_offset = 0
            plt.subplots_adjust(hspace=0)

            axs[0].scatter(
                df_all["eps1,eps2"] + constant_offset, df_all["phi_mean"], c=clr
            )
            axs[0].set_ylabel(r"$\phi$")
            axs[0].set_ylim(0.5, 1.0)
            axs[0].set_yticks([0.6, 0.8, 1.0])

            axs[1].scatter(
                df_all["eps1,eps2"] + constant_offset, df_all["shape_mean"], c=clr
            )
            axs[1].set_ylabel(r"$C$")
            axs[1].set_ylim(0.77, 0.93)
            axs[1].set_yticks([0.8, 0.84, 0.88, 0.92])

            axs[2].errorbar(
                df_all["eps1,eps2"] + df_all["Offset"] + constant_offset,
                df_all["v_mean"],
                yerr=df_all["v_std"],
                fmt="o",
                capsize=5,
                c=clr,
            )
            axs[2].set_ylabel(r"v $(\mu m/min)$")
            axs[2].set_ylim(0, 1)

            axs[3].scatter(
                df_all["eps1,eps2"] + constant_offset, df_all["NE_mean"], c=clr
            )
            axs[3].set_ylabel(r"$NE (cell \cdot min)^{-1}$")
            axs[3].set_ylim(0, 0.22)

            axs[3].set_xticks(range(len(df_all["params"].unique())))
            xticklabels = [
                full_label.split()[0] + " " + full_label.split()[1][:-1]
                for full_label in df_all["params"].unique()
            ]
            axs[3].set_xticklabels(xticklabels)
            plt.tight_layout()
            plt.savefig(
                f"simulationStackedObservablesGroupedk_off{fixed_val}_kl={kl}.png"
            )
            # plt.close(fig)

            debugpy.breakpoint()

            # plt.show()

        # plt.show()


if __name__ == "__main__":
    main()
