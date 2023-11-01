import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import debugpy
import csv
from io import StringIO

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
        # omit last row, which represents contacts with the boundary
        matrix = matrix[:-1, :-1]

        matrices.append(matrix)
    return matrices


def calculateNeighborExchanges(contactMatrix):
    diff = np.zeros(len(contactMatrix)-1)
    for i in range(len(contactMatrix)-1):
        diff[i] = np.count_nonzero(contactMatrix[i+1] - contactMatrix[i])
        nonzeros = np.nonzero(contactMatrix[i+1] - contactMatrix[i])
    return diff

def process_data(att, v0, att2, seed, fileheader):
    # fileheader = "test6"
    outputFileheader = "output" + fileheader[8:]
    fileExtensions = [".xStream", ".xMinStream",
                      ".cijStream", ".cijMinStream", ".comStream", ".shapeStream"]
    filenames = [fileheader + ext for ext in fileExtensions]
    # coordinates of all vertices
    cellPos = np.loadtxt(filenames[0])
    # energy minimized coordinates
    #cellMinPos = np.loadtxt(filenames[1])
    #cij = readCSV(filenames[2])  # cell-cell contact network
    # contact of minimized coordinates
    #cijMin = readCSV(filenames[3])
    #cellCOM = np.loadtxt(filenames[4])
    cellShape = np.loadtxt(filenames[5])

    # extract header information from cellPos and then remove header
    #NCELLS = int(cellPos[0][0])
    #NVTOT = int(cellPos[0][1])
    cellPos = cellPos[1:]
    #globalXLim = [min(cellPos[:, 0]), max(cellPos[:, 0])]
    #globalYLim = [min(cellPos[:, 1]), max(cellPos[:, 1])]
    shapeParameters = np.ravel(cellShape[:, 0::2])

    data_shape = {
        'shapeParameter': shapeParameters,
        'att': float(att),
        'v0': float(v0),
        'att2': float(att2),
        'seed': seed
    }
    df_shape = pd.DataFrame(data_shape)

    #minNE = calculateNeighborExchanges(cijMin)
    minNE = 0
    data_NE = {
        'minNE': minNE,
        'att': float(att),
        'v0': float(v0),
        'att2': float(att2),
        'seed': seed
    }
    #df_NE = pd.DataFrame(data_NE)

    #return df_shape, df_NE
    return df_shape


def main():
    att_arr = ["0.0006", "0.006", "0.06"]
    v0_arr = ["0.1"]
    att2_arr = ["0.0012", "0.012"]
    t_abp = "1.0"
    phi = "0.8"
    seeds = 10

    # psm_calA01.15_phi0.8_tm0_v00.1_t_abp1.0/_N40_dur200_att0.06_att20.012_start1_end10_sd10

    df_shapes = pd.DataFrame() 
    df_NEs = pd.DataFrame()
    # load packing fraction data into dataframe
    df_phis = pd.read_csv(f"C:\\Users\\atata\\projects\\psm_extracellular_calculation\\windowedPhiDataFrame.txt")

    # load shape and neighbor exchange data into dataframes
    for att in att_arr:
        for v0 in v0_arr:
            for att2 in att2_arr:
                for seed in range(1, seeds+1):
                    fileheader = f"pipeline\\cells\\psm\\psm_calA01.15_phi{phi}_tm0_v0{v0}_t_abp{t_abp}\\_N40_dur200_att{att}_att2{att2}_start1_end{seeds}_sd{seed}"
                    #df_shape, df_NE = process_data(att, v0, att2, seed, fileheader)
                    df_shape = process_data(att, v0, att2, seed, fileheader)
                    df_shapes = pd.concat([df_shapes, df_shape])
                    #df_NEs = pd.concat([df_NEs, df_NE])

    # Perform grouping on dataframes in order to make summary plots.

    # Add a column representing unique combinations of parameters. 
    # Use lambda here to convert combo column values into strings,
    #  then join the strings into one string so it can be read through seaborn
    # group the dataframe by combo parameters, then average over matching seeds and return a dataframe

    #df_NEs["(att, att2, v0)"] = df_NEs.apply(lambda row: ', '.join(
    #    [str(row['att']), str(row['att2']), str(row['v0'])]), axis=1)    
    #df_NEs_means = df_NEs.groupby(["(att, att2, v0)", "seed"]).mean().reset_index()

    df_shapes["(att, att2, v0)"] = df_shapes.apply(lambda row: ', '.join(
        [str(row['att']), str(row['att2']), str(row['v0'])]), axis=1)
    df_shapes_means = df_shapes.groupby(
        ["(att, att2, v0)", "seed"]).mean().reset_index()
    
    df_phis["(att, att2, v0)"] = df_phis.apply(lambda row: ', '.join(
        [str(row['att']), str(row['att2']), str(row['v0'])]), axis=1)
    df_phis_means = df_phis.groupby(
        ["(att, att2, v0)", "seed"]).mean().reset_index()

    # make boxplots of the mean values for each seed for each parameter combo
    #sns_NEs = sns.catplot(df_NEs_means, x="(att, att2, v0)", y="minNE", kind="box", color="skyblue", height = 8, aspect=2)
    sns_shapes = sns.catplot(df_shapes_means, x="(att, att2, v0)",
                y="shapeParameter", kind="box", color="skyblue", height=8, aspect=2)
    sns_phis = sns.catplot(df_phis_means, x="(att, att2, v0)",
                y="phi", kind="box", color="skyblue", height=8, aspect=2)
    
    #sns_NEs.figure.savefig(f"sns_NEs_{v0_arr[0]}.png")
    sns_shapes.figure.savefig(f"sns_shapes_{v0_arr[0]}.png")
    sns_phis.figure.savefig(f"sns_phis_{v0_arr[0]}.png")

    
    # make boxplots for each parameter combo, pooling observations from each seed
    #sns.catplot(df_NEs, x="(att, att2, v0)",
    #            y="minNE", kind="box", color="skyblue", height=8, aspect=2)
    sns.catplot(df_shapes, x="(att, att2, v0)",
                y="shapeParameter", kind="box", color="skyblue", height=8, aspect=2)
    sns.catplot(df_phis, x="(att, att2, v0)",
                y="phi", kind="box", color="skyblue", height=8, aspect=2)

    debugpy.breakpoint()

    plt.show()


if __name__ == "__main__":
    main()
