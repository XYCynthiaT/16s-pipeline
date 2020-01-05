import zipfile
import pandas as pd
import csv

def unzip(zipfn, destination):
    with zipfile.ZipFile(zipfn, 'r') as zip_ref:
        zip_ref.extractall(destination)
    
zipfile1 = snakemake.input[0]
zipfile2 = snakemake.input[1]

unzip(zipfile1, "output/s4FastQC")
unzip(zipfile2, "output/s4FastQC")

paths = ["output/s4FastQC/r1_fastqc/fastqc_data.txt", "output/s4FastQC/r2_fastqc/fastqc_data.txt"]

def readtbl(path):
    '''
    select rows by modules
    '''
    tbl = pd.read_csv(path, sep = "\t", header = 12, nrows = 57)
    return tbl

# def cleanColnames(tbl):
#     tbl.columns = tbl.columns.str.strip().str.replace(' ', '_')
#     return tbl

def cutLength(tbl):
    cutLen = 200
    for i in range(len(tbl.index)):
        if float(tbl.iloc[i][5]) <= 30:
            positions = tbl.iloc[i][0]
            cutLen = positions.split("-")[0]
            break
    print(cutLen)
    return int(cutLen)
            
cutLens = []
for path in paths:
    scores = readtbl(path)
    cutLens.append(cutLength(scores))

with open('cutLens.csv', 'w') as myfile:
        writer = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        writer.writerow(['r1', 'r2'])
        writer.writerow(cutLens)

