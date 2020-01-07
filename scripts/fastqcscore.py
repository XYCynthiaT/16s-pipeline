import zipfile
import pandas as pd
import csv
import io
from snakemake import shell


def unzip(zipfn, destination):
    with zipfile.ZipFile(zipfn, 'r') as zip_ref:
        zip_ref.extractall(destination)

def readtbl(path):
    #FIXME select rows by modules
    # tbl = pd.read_csv(path, sep = "\t", header = 12, nrows = 57)
    # return tbl
    path = "output/s4FastQC/r1_fastqc/fastqc_data.txt"
    with open(path, "rt") as fh:
        start = False
        tbl = ""
        for l in fh:
            if start:
                if l.startswith(">>END_MODULE"):
                    break
                tbl += l
                continue
            if l.startswith(">>Per base sequence quality"):
                start = True
    return pd.read_csv(io.StringIO(tbl), sep = "\t")

def cutLength(tbl):
    cutLen = 200
    for i in range(len(tbl.index)):
        if float(tbl.iloc[i][5]) <= 30:
            positions = tbl.iloc[i][0]
            cutLen = positions.split("-")[0]
            break
    print(cutLen)
    return int(cutLen)

def main():
    zipfile1 = snakemake.input[0]
    zipfile2 = snakemake.input[1]

    unzip(zipfile1, "output/s4FastQC")
    unzip(zipfile2, "output/s4FastQC")

    paths = ["output/s4FastQC/r1_fastqc/fastqc_data.txt", 
             "output/s4FastQC/r2_fastqc/fastqc_data.txt"]

    cutLens = []
    
    for path in paths:
        scores = readtbl(path)
        cutLens.append(cutLength(scores))

    with open('cutLens.csv', 'w') as myfile:
            writer = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            writer.writerow(['r1', 'r2'])
            writer.writerow(cutLens)
            

if __name__ == "__main__":
    main()