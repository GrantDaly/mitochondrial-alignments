import pdb
import pandas as pd
import numpy as np
from pathlib import Path
import pyranges as pr
import sys
from pathlib import Path
import json
import glob
import datetime
# insert stats

#

def getRawCoverageDF(bedPrefix, inputDir, aggregate=False):
    covGlob = glob.glob(str(inputDir) + "/*." + bedPrefix + ".coverage.tsv")
    #pdb.set_trace()
    covDtypes = {
    "Chromosome":"object",
    "Start":"int64",
    "End":"int64",
    "Name":"object",
    "Score":"float64",
    "Strand":"object",
    "Sample":"object",
    "Offset":"int64",
    "Forward_Depth": "int64",
    "Reverse_Depth": "int64",
        "Depth":"int64"}

    coverageList = []
    
    for filename in covGlob:
        tempCovDF = pd.read_csv(filename, sep="\t", dtype=covDtypes)
        if(aggregate == True):
            
            tempCovDF = tempCovDF.groupby(["Name", "Chromosome", "Start", "End", "Strand","Score","Sample"])["Depth"].agg('sum').reset_index()
            tempCovDF = tempCovDF.rename(columns={'Depth': 'Total Depth'})
            tempCovDF['Length'] = tempCovDF['End'] - tempCovDF['Start']
            tempCovDF['Mean Depth'] = tempCovDF['Total Depth'] / tempCovDF['Length']
        coverageList.append(tempCovDF)


    covGlobStats = glob.glob(str(inputDir) + "/*." + bedPrefix + ".coverage.stats.tsv")
    coverageStatList = []
    for statsFilename in covGlobStats:
        tempCovStatsDF = pd.read_csv(statsFilename, sep="\t")
        coverageStatList.append(tempCovStatsDF)

    outCovDF = pd.concat(coverageList, axis=0)
    outCovStatsDF = pd.concat(coverageStatList, axis=0)

    return outCovDF, outCovStatsDF

def getRawInsertsDF(bedPrefix, inputDir):

    insertDtypes = {"Sample":"object",
                "Intervals": "int64",
	        "Raw Density":"float64", "% Density":"float64",
	        "Smoothed % Density":"float64"}

    insertHistList = []
    insertHistGlob = glob.glob(str(inputDir) + "/*." + bedPrefix + ".insert.hist.tsv")
    for insertFileName in insertHistGlob:
        tempHist = pd.read_csv(insertFileName, sep="\t", dtype=insertDtypes)
        insertHistList.append(tempHist)
    insertHistDF = pd.concat(insertHistList, axis=0)

    insertStatsList = []
    insertStatsGlob = glob.glob(str(inputDir) + "/*." + bedPrefix + ".insert.stats.tsv")
    for insertStatFile in insertStatsGlob:
        tempStats = pd.read_csv(insertStatFile, sep="\t")
        insertStatsList.append(tempStats)
    insertStatDF = pd.concat(insertStatsList, axis=0)

    return insertHistDF, insertStatDF
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='generate analysis workbooks from raw files')

    parser.add_argument('-p', '--parameters')
    parser.add_argument('-m', '--metadata')
    parser.add_argument('-i', '--inputdir')
    parser.add_argument('-o', '--outputdir')
    parser.add_argument('-n', '--nameprefix', default="DAMPs")
    parser.add_argument('-t', '--timestamp', action='store_true')

    args = parser.parse_args()
    #print(args)

    with open(args.parameters) as paramFile:
        params = json.load(paramFile)

    with open(args.metadata) as metaFile:
        metadata = pd.read_csv(metaFile, sep="\t")

    inputDir = Path(args.inputdir)
    outputDir = Path(args.outputdir)

    # go through groups and read in all requested files
    coverageDict = {}
    coverageStatsDict = {}
    insertHistDict = {}
    insertStatsDict = {}
    
    for bedParams in params['groups']:

        tempName = bedParams['name']
        
        if(bedParams['positions'] == "full"):
            # don't need to aggregate
            tempCov , tempCovStats = getRawCoverageDF(tempName, inputDir, aggregate=False)
            coverageDict[tempName] = tempCov
            coverageStatsDict[tempName] = tempCovStats
            
        elif(bedParams['positions'] == "aggregate"):
            # aggregate by interval
            tempCov , tempCovStats = getRawCoverageDF(tempName, inputDir, aggregate=True)
            coverageDict[tempName] = tempCov
            coverageStatsDict[tempName] = tempCovStats

            
        # if inserts are requested, add to dict
        if(bedParams['inserts'] == "include"):
            # add to dict
            tempInsertHist, tempInsertStats = getRawInsertsDF(tempName, inputDir)
            insertHistDict[tempName] = tempInsertHist
            insertStatsDict[tempName] = tempInsertStats
            
        elif(('inserts' not in bedParams.keys()) or (bedParams['inserts'] == "exclude")):
            pass

    # after the raw coverages and inserts are read in, find the mitochondrial value and numt value, and normalize them.
        # will normalize the 1bp resolution set and the mean of the stats set
    mitoName = params['normalization']['mitochondrial-name']
    numtName = params['normalization']['numt-name']
    
    # start with normalizing stats df, for this I'm only doing mitochondrial and NUMT
    normMitoStatsDF = coverageStatsDict[mitoName]

    # use the mean NUMT depth to normalize coverage. Will repeat this for the 1bp or interval-level coverages as well
    normNumtDepth = coverageStatsDict[numtName].loc[:, ["Sample", "Mean"]]
    normNumtDepth = normNumtDepth.rename(columns={"Mean": "NUMT Mean"})

    normMitoStatsDF = normMitoStatsDF.merge(normNumtDepth, on="Sample")
    normMitoStatsDF['Norm Mean'] = normMitoStatsDF['Mean'] / normMitoStatsDF['NUMT Mean'] 
    normMitoStatsDF = normMitoStatsDF.drop(columns="NUMT Mean")

    # replace the mitochondrial df with the updated one
    coverageStatsDict.pop(mitoName)
    coverageStatsDict[mitoName] = normMitoStatsDF
    
    # for mito, NUMT, Shuffled, etc. normalize by the mean NUMT for that sample.
    for dfName, covDF in coverageDict.items():
        
        covDF = covDF.merge(normNumtDepth, how="left", on="Sample")
        if("Depth" in covDF.columns):
            covDF['Norm Depth'] = covDF['Depth'] / covDF['NUMT Mean']
        elif("Mean Depth" in covDF.columns):
            covDF['Norm Depth'] = covDF['Mean Depth'] / covDF['NUMT Mean']

        # once I'm confident in results could drop the column
        #covDF = covDF.drop(columns="NUMT Mean")
        coverageDict[dfName] = covDF

    #### for mitochondria, add a 100 bp binned sheet ######
    binnedMito = coverageDict[mitoName].loc[:,["Sample", "Offset","Depth", "Norm Depth"]]
    binnedMito['Bin'] = (binnedMito['Offset'] // 100 ) + 1
    binnedMito = pd.pivot_table(binnedMito, values=["Depth", "Norm Depth"], index=["Bin"], columns="Sample", aggfunc="mean")

    # process inserts
    # insert stats should be ready

    # insert histograms need to be pivoted
    pivotedInsertHistDict = {}
    for insertName, insertHist in insertHistDict.items():
        tempPivot = insertHist.pivot_table(index=["Intervals"], columns=["Sample"],
                                                  values=["Smoothed % Density"])
        pivotedInsertHistDict[insertName] = tempPivot

    # start writing excel Workbook
    if args.timestamp == True:
        now = datetime.datetime.now()
        outExcelName = args.nameprefix + "-" +  str(now) + ".xlsx"
    else:
        outExcelName = args.nameprefix + ".xlsx"

    outDirWorkbooks = Path(outputDir  / "workbooks")
    if not outDirWorkbooks.is_dir():
        outDirWorkbooks.mkdir()

    writer = pd.ExcelWriter(outDirWorkbooks / outExcelName)

    for tempName, tempDF in coverageStatsDict.items():
        tempDF.to_excel(writer, tempName + " Coverage Stats",index=None)

    for tempName, tempDF in coverageDict.items():
        tempDF.to_excel(writer, tempName + " Coverage",index=None)

    binnedMito.to_excel(writer, "100bp Bin Mito")
    
    for tempName, tempDF in insertStatsDict.items():
        tempDF.to_excel(writer, tempName + " Insert Stats",index=None)
    for tempName, tempDF in pivotedInsertHistDict.items():
        tempDF.to_excel(writer, tempName + " Insert Histograms")

    writer.save()
