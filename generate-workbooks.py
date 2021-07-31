import pdb
import pandas as pd
import numpy as np
from pathlib import Path
import pyranges as pr
import sys
from pathlib import Path
# insert stats

workingDir = Path(sys.argv[1])
designDF = pd.read_csv(workingDir / "design.csv", sep="\t")

insertDtypes = {"Chromosome":"object",
                "Start": "int64","End":"int64",
	        "Name":"object","Score":"float64","Strand":"object",
	        "Sample":"object","Offset":"int64","Depth":"int64"}
insertDir = Path(workingDir / "insert-outputs")
insertDtypes = {}

numtInsertStats = pd.read_csv(insertDir / "Numt-InsertStats.tsv", sep="\t",dtype=insertDtypes)

numtInsertHist = pd.read_csv(insertDir / "Numt-InsertHist.tsv", sep="\t",dtype=insertDtypes)


mitoInsertStats= pd.read_csv(insertDir / "Mito-InsertStats.tsv", sep="\t",dtype=insertDtypes)

mitoInsertHist= pd.read_csv(insertDir / "Mito-InsertHist.tsv", sep="\t",dtype=insertDtypes)


numtInsertStats["Insert Origin"] = "NUMT"
numtInsertHist["Insert Origin"] = "NUMT"
mitoInsertStats["Insert Origin"] = "Mito"
mitoInsertHist["Insert Origin"] = "Mito"

insertStatsAll = pd.concat([numtInsertStats, mitoInsertStats], axis=0)
insertHistAll = pd.concat([numtInsertHist, mitoInsertHist], axis=0)

insertStatsAll = insertStatsAll.merge(designDF, on="Sample", how="left")
insertHistAll = insertHistAll.merge(designDF, on="Sample", how="left")

covDtypes = {
    "Chromosome":"object",
    "Start":"int64",
    "End":"int64",
    "Name":"object",
    "Score":"float64",
    "Strand":"object",
    "Sample":"object",
    "Offset":"int64",
    "Depth":"int64"
}
# coverage stats
coverageDir = Path(workingDir / "coverage-outputs")

mitoCovCutdown = pd.read_csv(coverageDir / "Mito-Coverages.tsv", sep="\t",dtype=covDtypes)

numtCovCutdown = pd.read_csv(coverageDir / "Numt-Coverages.tsv", sep="\t",dtype=covDtypes)

shufNumtCovCutdown = pd.read_csv(coverageDir / "Shuffled-Coverages.tsv", sep="\t",dtype=covDtypes)

dayamaCovCutdown = pd.read_csv(coverageDir / "Dayama-Coverages.tsv", sep="\t",dtype=covDtypes)
# these are the coverages aligned only to the mitochondria


#aggregate by numt
covByNumt = numtCovCutdown.groupby(["Name", "Chromosome", "Start", "End", "Strand","Score","Sample"])["Depth"].agg('sum').reset_index()
covByNumt['Length'] = covByNumt['End'] - covByNumt['Start']
covByNumt['Mean'] = covByNumt['Depth'] / covByNumt['Length']

#aggregate by shuffle
covByShuf = shufNumtCovCutdown.groupby(["Name", "Chromosome", "Start", "End", "Strand","Score","Sample"])["Depth"].agg('sum').reset_index()
covByShuf['Length'] = covByShuf['End'] - covByShuf['Start']
covByShuf['Mean'] = covByShuf['Depth'] / covByShuf['Length']

#aggregate by Dayama
covByDayama = dayamaCovCutdown.groupby(["Name", "Chromosome", "Start", "End", "Strand","Score","Sample"])["Depth"].agg('sum').reset_index()
covByDayama['Length'] = covByDayama['End'] - covByDayama['Start']
covByDayama['Mean'] = covByDayama['Depth'] / covByDayama['Length']


# normalize each site in mito by dividing by mean NUMT coverage
numtStats = numtCovCutdown.loc[:, ["Sample", "Depth"]].groupby("Sample").describe()
numtMeans = numtCovCutdown.loc[:, ["Sample", "Depth"]].groupby("Sample").agg('mean')

# Normalize by mean numt coverages

humanMtNorm = mitoCovCutdown.merge(numtMeans, how="left", on="Sample", suffixes=["_Mito", "_NUMT"])
humanMtNorm["Norm Depth"] = humanMtNorm["Depth_Mito"] / humanMtNorm["Depth_NUMT"]


covByNumt = covByNumt.merge(numtMeans.rename(columns={"Depth": "Mean All NUMT"})
                    , how="left", on="Sample")
covByNumt['Normalized Mean'] = covByNumt['Mean'] / covByNumt['Mean All NUMT']

covByShuf = covByShuf.merge(numtMeans.rename(columns={"Depth": "Mean All NUMT"})
                    , how="left", on="Sample")
covByShuf['Normalized Mean'] = covByShuf['Mean'] / covByShuf['Mean All NUMT']

covByDayama = covByDayama.merge(numtMeans.rename(columns={"Depth": "Mean All NUMT"})
                    , how="left", on="Sample")
covByDayama['Normalized Mean'] = covByDayama['Mean'] / covByDayama['Mean All NUMT']

# bin 100 bp bins
humanMtNormBins = humanMtNorm.reset_index().copy()
humanMtNormBins['Bin'] = (humanMtNormBins['Offset'] // 100 ) + 1
binnedMeans = pd.pivot_table(humanMtNormBins, values=["Depth_Mito", "Norm Depth"], index=["Bin"], columns="Sample", aggfunc="mean")

# mitochondrial stats
mitoStats = humanMtNorm.groupby('Sample')[['Depth_Mito', 'Norm Depth', 'Depth_NUMT']].describe()
mitoStatsFlat = mitoStats.copy()
mitoStatsFlat.columns = mitoStatsFlat.columns.map('_'.join)

# mitoNormStats = humanMtNorm.groupby('Sample')['Norm Depth'].describe()
mitoStatsFlat['Depth_Mito_cv'] = mitoStatsFlat['Depth_Mito_mean'] / mitoStatsFlat['Depth_Mito_std']
mitoStatsFlat['Norm Depth_cv'] = mitoStatsFlat['Norm Depth_mean'] / mitoStatsFlat['Norm Depth_std']
mitoStatsFlat['Depth_NUMT_cv'] = mitoStatsFlat['Depth_NUMT_mean'] / mitoStatsFlat['Depth_NUMT_std']

meanCovMerged = mitoStatsFlat.merge(designDF, how="left", on="Sample")

# write out to workbook
outDirWorkbooks = Path(workingDir / "workbooks")
if not outDirWorkbooks.is_dir():
    print("make directory")
    outDirWorkbooks.mkdir()
import datetime
outName = sys.argv[2]
now = datetime.datetime.now()
writer = pd.ExcelWriter(outDirWorkbooks / ('DAMP-Stats-' + outName + str(now) + '.xlsx'), engine='xlsxwriter')
insertStatsAll.to_excel(writer, "Insert Stats",index=False)

insertHistPivot = insertHistAll.pivot_table(index=["Intervals"], columns=["Sample","Insert Origin"], values=["% Density"])
insertHistPivot.to_excel(writer, "Insert Distributions")

insertHistPivot = insertHistAll.pivot_table(index=["Intervals"], columns=["Sample","Insert Origin"], values=["Smoothed % Density"])
insertHistPivot.to_excel(writer, "Smoothed Insert Distributions")

meanCovMerged.to_excel(writer, "Coverage Stats",index=False)

#binned mito coverage
binnedMeans.to_excel(writer, "100bp Binned Coverage")

#pivot mito coverage
covPivot = humanMtNorm.pivot_table(index=["Chromosome", "Start", "End", "Name", "Score", "Strand", "Offset"],
                 columns="Sample", values=["Depth_Mito", "Depth_NUMT", "Norm Depth"]).reset_index()
covPivot.to_excel(writer, "Mito Coverage")

#pivot numt
covNumtPivot = covByNumt.pivot_table(index=["Chromosome", "Start", "End", "Name", "Score", "Strand", "Length"],
                 columns="Sample", values=["Depth","Length", "Mean","Mean All NUMT","Normalized Mean"]).reset_index()
covNumtPivot.to_excel(writer, "Numt Mean Coverage")

covShufPivot = covByShuf.pivot_table(index=["Chromosome", "Start", "End", "Name", "Score", "Strand", "Length"],
                 columns="Sample", values=["Depth","Length", "Mean","Mean All NUMT","Normalized Mean"]).reset_index()
covShufPivot.to_excel(writer, "Shuf Mean Coverage")

covDiyamaPivot = covByDayama.pivot_table(index=["Chromosome", "Start", "End", "Name", "Score", "Strand", "Length"],
                 columns="Sample", values=["Depth","Length", "Mean","Mean All NUMT","Normalized Mean"]).reset_index()
covDiyamaPivot.to_excel(writer, "Dayama Mean Coverage")

writer.save()
 
# generate bedgraphs and bigwigs
outDirRaw = Path(coverageDir / "raw-coverages")
outDirNorm = Path(coverageDir / "norm-coverages")
if not outDirRaw.is_dir():
    outDirRaw.mkdir()
if not outDirNorm.is_dir():
    outDirNorm.mkdir()
    "Human-Numt-Coverage-Raw.tsv"

for name, group in humanMtNorm.groupby("Sample"):
    outGroup = group.copy().reset_index(drop=True)
    outGroup['Start'] = outGroup['Offset']
    outGroup['End']  = outGroup['Start'] + 1
    outGroupNorm = outGroup.copy()
    outGroup['Count'] = outGroup['Depth_Mito']
    outGroupNorm['Count'] = outGroup['Norm Depth']

    outNameRaw = name + "-raw"
    outNameNorm = name + "-norm"

#     filenameList.append(outName)

    outGroup = outGroup.loc[:, ["Chromosome", "Start", "End", "Count"]]
    outGroup.to_csv(outDirRaw / (outNameRaw + ".bg"), sep="\t", index=False, header=False)
    
    outGroupNorm = outGroupNorm.loc[:, ["Chromosome", "Start", "End", "Count"]]
    outGroupNorm.to_csv(outDirNorm / (outNameNorm + ".bg"), sep="\t", index=False, header=False)

    # pyranges to_bigwig not working
    # tempPyRanges = pr.PyRanges(outGroup)
    # #pdb.set_trace()
    # outRawString = str(outDirRaw / (outNameRaw + ".bw"))
    # tempPyRanges.to_bigwig(path = outRawString)
    
    # tempPyRangesNorm = pr.PyRanges(outGroupNorm)
    # outNormString = str(outDirNorm / (outNameNorm + ".bw"))
    # tempPyRangesNorm.to_bigwig(path = outNormString)
