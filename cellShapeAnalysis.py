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


def main():
    att_arr = ["0.001", "0.01", "0.1"]
    #att_arr = ["0.001", "0.01"]
    v0_arr = ["0.02"]
    #v0_arr = ["0.02", "0.04"]
    #k_ecm_arr = ["0.05", "0.5"]
    k_ecm_arr = ["0.05", "0.5", "5"]
    seeds = 10
    df_shapes = pd.DataFrame()
    df_NEs = pd.DataFrame()
    for att in att_arr:
        for v0 in v0_arr:
            for k_ecm in k_ecm_arr:
                for seed in range(1, seeds+1):
                    # fileheader = "test6"
                    fileheader = "pipeline\cells\psm\psm_calA01.0_phi0.74_tm10.0_v0"+v0+"_t_abp1.0k_ecm"+k_ecm+"k_off1.0\_N40_dur250_att"+att+"_start1_end"+str(seeds)+"_sd"+str(seed)
                    outputFileheader = "output"+fileheader[8:]
                    fileExtensions = [".xStream", ".xMinStream",
                                    ".cijStream", ".cijMinStream", ".comStream", ".shapeStream"]
                    filenames = [fileheader + ext for ext in fileExtensions]
                    cellPos = np.loadtxt(filenames[0])  # coordinates of all vertices
                    cellMinPos = np.loadtxt(filenames[1])  # energy minimized coordinates
                    cij = readCSV(filenames[2])  # cell-cell contact network
                    # cell-cell contact network of minimized coordinates
                    cijMin = readCSV(filenames[3])
                    cellCOM = np.loadtxt(filenames[4])
                    cellShape = np.loadtxt(filenames[5])

                    # extract header information from cellPos and then remove header
                    NCELLS = int(cellPos[0][0])
                    NVTOT = int(cellPos[0][1])
                    cellPos = cellPos[1:]
                    globalXLim = [min(cellPos[:, 0]), max(cellPos[:, 0])]
                    globalYLim = [min(cellPos[:, 1]), max(cellPos[:, 1])]
                    shapeParameters = np.ravel(cellShape[:, 0::2])
                    df_shape = pd.DataFrame(shapeParameters, columns = ['shapeParameter'])
                    df_shape['att'] = float(att)
                    df_shape['v0'] = float(v0)
                    df_shape['k_ecm'] = float(k_ecm)
                    df_shape['seed'] = seed
                    df_shapes = pd.concat([df_shapes, df_shape])

                    minNE = calculateNeighborExchanges(cijMin)
                    df_NE = pd.DataFrame(minNE, columns = ['minNE'])
                    df_NE['att'] = float(att)
                    df_NE['v0'] = float(v0)
                    df_NE['k_ecm'] = float(k_ecm)
                    df_NE['seed'] = seed
                    df_NEs = pd.concat([df_NEs, df_NE])

    # add a column representing unique combinations of parameters. 
    #   use lambda here to convert combo column values into strings,
    #       then join the strings into one string so it can be read through seaborn
    df_NEs["(att, k_ecm, v0)"] = df_NEs.apply(lambda row: ', '.join(
        [str(row['att']), str(row['k_ecm']), str(row['v0'])]), axis=1)
    
    # group the dataframe by combo parameters, then average over matching seeds and return a dataframe
    df_NEs_means = df_NEs.groupby(["(att, k_ecm, v0)", "seed"]).mean().reset_index()

    df_shapes["(att, k_ecm, v0)"] = df_shapes.apply(lambda row: ', '.join(
        [str(row['att']), str(row['k_ecm']), str(row['v0'])]), axis=1)
    df_shapes_means = df_shapes.groupby(
        ["(att, k_ecm, v0)", "seed"]).mean().reset_index()

    plt.figure()
    sns.catplot(df_NEs_means, x="(att, k_ecm, v0)", y="minNE", kind="box", color="skyblue")

    plt.figure()
    sns.catplot(df_shapes_means, x="(att, k_ecm, v0)",
                y="shapeParameter", kind="box", color="skyblue")




if __name__ == "__main__":
    main()
