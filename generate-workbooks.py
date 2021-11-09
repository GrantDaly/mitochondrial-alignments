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
import gzip
import requests

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
    outCovStatsDF = pd.concat(coverageStatList, axis=0).sort_values("Sample")

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
    insertStatDF = pd.concat(insertStatsList, axis=0).sort_values("Sample")

    return insertHistDF, insertStatDF

def createBigWigs(coverageDir, designDF, normMitoCov):
    
    # generate bedgraphs and bigwigs
    outDirRaw = Path(coverageDir / "raw-coverages")
    outDirNorm = Path(coverageDir / "norm-coverages")
    if not outDirRaw.is_dir():
        outDirRaw.mkdir()
    if not outDirNorm.is_dir():
        outDirNorm.mkdir()
        "Human-Numt-Coverage-Raw.tsv"

    rawDesignList = []
    normDesignList = []
    for name, group in normMitoCov.groupby("Sample"):

        designEntry = designDF.query('Sample == @name')
        individual = designEntry['Individual ID'].to_list()[0]
        sampleName = designEntry['Sample'].to_list()[0]
        condition = designEntry['Condition'].to_list()[0]

        outGroup = group.copy().reset_index(drop=True)
        outGroup['Start'] = outGroup['Offset']
        outGroup['End']  = outGroup['Start'] + 1
        outGroupNorm = outGroup.copy()

    #    outGroup['Count'] = outGroup['Depth'] 
    #    outGroupNorm['Count'] = outGroup['Norm Depth']

        outNameRaw = name + "-raw"
        outNameNorm = name + "-norm"

    #     filenameList.append(outName)

        # only outputting bigwigs for now.

        outGroupRawAll = outGroup.loc[:, ["Chromosome", "Start", "End", "Depth", "Forward_Depth", "Reverse_Depth"]]

        #outGroupRawAll['Percent_Diff_Strand'] = ((outGroupRawAll['Forward_Depth'] - outGroupRawAll['Reverse_Depth']) /  outGroupRawAll['Reverse_Depth']) * 100
        #pdb.set_trace()
        outGroupRawAll['LogFC_Diff_Strand'] = np.log2((outGroupRawAll['Forward_Depth'] + 1) /  (outGroupRawAll['Reverse_Depth'] + 1))

        #outGroup.to_csv(outDirRaw / (outNameRaw + ".bg"), sep="\t", index=False, header=False)

        outGroupNorm = outGroupNorm.loc[:, ["Chromosome", "Start", "End", "Norm Depth"]]
        #outGroupNorm.to_csv(outDirNorm / (outNameNorm + ".bg"), sep="\t", index=False, header=False)

        # pyranges to_bigwig not working
        tempPyRangesRaw = pr.PyRanges(outGroupRawAll)

        outRawString = str(outDirRaw / (outNameRaw + ".bw"))
        outRawStringFor = str(outDirRaw / (outNameRaw + ".F1R2.bw"))
        outRawStringRev = str(outDirRaw / (outNameRaw + ".F2R1.bw"))
        outRawStringDiff = str(outDirRaw / (outNameRaw + ".StrandDiff.bw"))
        tempPyRangesRaw.to_bigwig(path = outRawString, value_col="Depth")
        #pdb.set_trace()
        rawDesignList.append({'TRACK_ID': outNameRaw + ".bw", 'INDIVIDUAL_ID': individual,
                              'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand': "Both",
                              'Condition-Strand': 'NA'})
        tempPyRangesRaw.to_bigwig(path = outRawStringFor, value_col="Forward_Depth")
        rawDesignList.append({'TRACK_ID': outNameRaw + ".F1R2.bw", 'INDIVIDUAL_ID': individual,
                              'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand':"Forward",
                              "Condition-Strand": condition + "-Forward"})

        tempPyRangesRaw.to_bigwig(path = outRawStringRev, value_col="Reverse_Depth")
        rawDesignList.append({'TRACK_ID': outNameRaw + ".F2R1.bw", 'INDIVIDUAL_ID': individual,
                              'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand':"Reverse",
                              'Condition-Strand': condition + "-Reverse"})

        tempPyRangesRaw.to_bigwig(path = outRawStringDiff, value_col="LogFC_Diff_Strand")
        rawDesignList.append({'TRACK_ID': outNameRaw + ".StrandDiff.bw", 'INDIVIDUAL_ID': individual,
                              'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand':"LogFC",
                              'Condition-Strand': "NA"})


        tempPyRangesNorm = pr.PyRanges(outGroupNorm)
        outNormString = str(outDirNorm / (outNameNorm + ".bw"))
        tempPyRangesNorm.to_bigwig(path = outNormString, value_col="Norm Depth")
        normDesignList.append({'TRACK_ID': outNameNorm + ".bw", 'INDIVIDUAL_ID': individual,
                              'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand':"Both"})


    outRawDesignDF = pd.DataFrame(rawDesignList)
    outRawDesignDF['DATA_Type'] = "Expression"

    outRawDesignDF.loc[:,["TRACK_ID", "DATA_Type", "INDIVIDUAL_ID", "SAMPLE_ID", "Condition","Strand","Condition-Strand"]]\
                             .to_csv(outDirRaw / "raw-attributes.txt", sep="\t", index=None)


    outNormDesignDF = pd.DataFrame(normDesignList)
    outNormDesignDF['DATA_Type'] = "Expression"

    outNormDesignDF.loc[:,["TRACK_ID", "DATA_Type", "INDIVIDUAL_ID", "SAMPLE_ID", "Condition","Strand"]]\
                             .to_csv(outDirNorm / "norm-attributes.txt", sep="\t", index=None)
    return

def makeRequest(idSeries):
    
        server = "https://rest.ensembl.org"
        ext = "/variation/homo_sapiens?genotypes=1&phenotypes=1"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        
        idList = idSeries.unique().tolist()
        idDict = {"ids": idList}
        r = requests.post(server+ext, headers=headers, data=json.dumps(idDict))

        if not r.ok:
          r.raise_for_status()
          sys.exit()

        decoded = r.json()
        return decoded

def getAnnotationsAll(idSeries):
    import requests, sys
 
    
#     r = requests.post(server+ext, headers=headers, data='{ "ids" : ["rs199476118", "rs879100564", "rs199476118" ] }')
    outList = []
    binSize=100
    idSeries = pd.Series(idSeries.unique())
    #pdb.set_trace()
    for rangeBin in range((len(idSeries) // binSize) + 1):
        decoded = makeRequest(idSeries[rangeBin*binSize:rangeBin*binSize + binSize])
        tempDF = processVEPReturn(decoded)
        outList.append(tempDF)
        
    decoded = makeRequest(idSeries[(rangeBin+1)*binSize:])
    tempDF = processVEPReturn(decoded)
    outList.append(tempDF)
    return pd.concat(outList,axis=0)


#     print(repr(decoded))
    return processVEPReturn(decoded)

def processVEPReturn(inJson):
    outList = []
    for rsid, entry in inJson.items():
#         print("###########" + rsid + "###########")
#         pprint(entry)
        entryDict = {}
        entryDict['ID'] = rsid
        if("clinical_significance" in entry.keys()):
            entryDict['clinical_significance'] = ":".join(entry['clinical_significance'])
        else:
            entryDict['clinical_significance'] = None
        
        if(len(entry['phenotypes']) > 0):
            phenotypeList = []
            for phenotypeEntry in entry['phenotypes']:
                phenotypeList.append(" ".join(["-".join((key,str(val))) for key,val in phenotypeEntry.items()]))
            entryDict['phenotypes'] = ":".join(phenotypeList)
        else:
            entryDict['phenotypes'] = None
    
        entryDict['most_severe_consequence'] = entry['most_severe_consequence']
        outList.append(entryDict)
    return pd.DataFrame(outList)


def parseVariantTSV(fileName):
    outList = []
    with gzip.open(fileName, 'rt') as inFile:
        for line in inFile:
            parsedLine = line.strip().split("\t")
            contig = parsedLine[0]
            locus = parsedLine[1]
            variantID = parsedLine[2]
            reference = parsedLine[3]
            alternate = parsedLine[4]
            varType = parsedLine[5]
            
            formatString = parsedLine[6]
            attributeNames = formatString.split(":")
            numAttributes = len(attributeNames)
            start = 7
            end = start + numAttributes
            
            numSamples = len(parsedLine[7:]) // numAttributes
            for _ in range(numSamples):
                sampleLine = parsedLine[start:end]
                start += numAttributes
                end += numAttributes
                
                entryDict = {"Chromosome": contig, "Start": locus, 
                             "ID": variantID,
                             "ref": reference, "alt": alternate, 
                                "Variant Type": varType}
                tempTuple = zip(attributeNames, sampleLine)
                for attr, val in tempTuple:
                    if(attr == "MinorFreq"):
                        refFreq, altFreq = val.split(",")
#                         entryDict["Ref AF"] = refFreq
                        entryDict["MinorFreq"] = altFreq
                    else:
                        entryDict[attr] = val
                outList.append(entryDict)
    outDF = pd.DataFrame(outList)
    
    outDF = outDF.replace('.', np.NaN)
    
    outDF = outDF.set_index(["Chromosome", "Start", "ref", "alt"])
    return outDF    



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='generate analysis workbooks from raw files')

    parser.add_argument('-p', '--parameters')
    parser.add_argument('-m', '--metadata')
    parser.add_argument('-i', '--inputdir')
    parser.add_argument('-o', '--outputdir')
    parser.add_argument('-n', '--nameprefix', default="DAMPs")
    parser.add_argument('-t', '--timestamp', action='store_true')
    parser.add_argument('-r', '--regenfiles', action='store_true')

    args = parser.parse_args()
    #print(args)

    with open(args.parameters) as paramFile:
        params = json.load(paramFile)

    #with open(args.metadata) as metaFile:
    #    metadata = pd.read_csv(metaFile, sep="\t")

    designDF = pd.read_csv(args.metadata, sep="\t")
    inputDir = Path(args.inputdir)
    outputDir = Path(args.outputdir)
    if not outputDir.is_dir():
        outputDir.mkdir()

    # go through groups and read in all requested files
    coverageDict = {}
    pivotCoverageDict = {}
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


    # replace the mitochondrial df with the updated one
    coverageStatsDict.pop(mitoName)
    coverageStatsDict[mitoName] = normMitoStatsDF
    
    # for mito, NUMT, Shuffled, etc. normalize by the mean NUMT for that sample.
    for dfName, rawCovDF in coverageDict.items():
        
        rawCovDF = rawCovDF.merge(normNumtDepth, how="left", on="Sample")
        if("Depth" in rawCovDF.columns):
            rawCovDF['Norm Depth'] = rawCovDF['Depth'] / rawCovDF['NUMT Mean']
            # for the 1bp resolution set, I just want to list the start position
            rawCovDF['Position'] = rawCovDF["Offset"] + 1
            
            collapsedCovDF = rawCovDF.pivot(index=["Position"],columns="Sample",values=['Forward_Depth', 'Reverse_Depth', 'Depth', 'Norm Depth'])
            
        elif("Mean Depth" in rawCovDF.columns):
            rawCovDF['Norm Depth'] = rawCovDF['Mean Depth'] / rawCovDF['NUMT Mean']
            collapsedCovDF = rawCovDF.pivot(index=["Chromosome", "Start", "End","Name"],columns="Sample",values=['Mean Depth', 'Norm Depth'])

        # once I'm confident in results could drop the column
        #covDF = covDF.drop(columns="NUMT Mean")
        coverageDict[dfName] = rawCovDF
        pivotCoverageDict[dfName] = collapsedCovDF
        


    #### for mitochondria, add a 100 bp binned sheet ######
    binnedMito = coverageDict[mitoName].loc[:,["Sample", "Offset","Depth", "Norm Depth"]]
    binnedMito['Bin'] = (binnedMito['Offset'] // 100 ) + 1
    binnedMito = pd.pivot_table(binnedMito, values=["Depth", "Norm Depth"], index=["Bin"], columns="Sample", aggfunc="mean")

    # process inserts
    # insert stats should be ready

    # insert histograms need to be pivoted
    
    outDirInserts = Path(outputDir  / "inserts")
    if not outDirInserts.is_dir():
        outDirInserts.mkdir()
    pivotedInsertHistDict = {}
    for insertName, insertHist in insertHistDict.items():
        insertHist.to_csv(outDirInserts / (insertName + ".inserts.csv"),sep="\t",index=None)
        tempPivot = insertHist.pivot_table(index=["Intervals"], columns=["Sample"],
                                                  values=["Smoothed % Density"])
        pivotedInsertHistDict[insertName] = tempPivot
    
    # start writing excel Workbook
    if args.timestamp == True:
        now = datetime.datetime.now().strftime("%m-%d-%Y-%H%M")
        outExcelName = args.nameprefix + "-" +  str(now) + ".xlsx"
        outExcelNameHeteroplasmy = args.nameprefix + "-Heteroplasmy-" +  str(now) + ".xlsx"
    else:
        outExcelName = args.nameprefix + ".xlsx"
        outExcelNameHeteroplasmy = args.nameprefix + "-Heteroplasmy.xlsx"

    outDirWorkbooks = Path(outputDir  / "workbooks")
    if not outDirWorkbooks.is_dir():
        outDirWorkbooks.mkdir()

    writer = pd.ExcelWriter(outDirWorkbooks / outExcelName)
    designDF.sort_values("Sample").to_excel(writer, "Design", index=None)

    for tempName, tempDF in coverageStatsDict.items():
        tempDF = tempDF.merge(designDF, how="left", on="Sample")
        tempDF.to_excel(writer, tempName + " Coverage Stats",index=None)

    for tempName, tempDF in pivotCoverageDict.items():
        tempDF.to_excel(writer, tempName + " Coverage")

    binnedMito.to_excel(writer, "100bp Bin Mito")
    
    for tempName, tempDF in insertStatsDict.items():
        tempDF = tempDF.merge(designDF, how="left", on="Sample")
        tempDF.to_excel(writer, tempName + " Insert Stats",index=None)
    for tempName, tempDF in pivotedInsertHistDict.items():
        tempDF.to_excel(writer, tempName + " Insert Histograms")

    writer.save()

    # write bigwig files
    outDirCoverage = Path(outputDir  / "coverage")
    if not outDirCoverage.is_dir():
        outDirCoverage.mkdir()
    # output raw coverage tsv
    for tempName, tempDF in coverageDict.items():
        tempDF.to_csv(outDirCoverage / (tempName + ".coverage.tsv"),sep="\t",index=None)
    
    createBigWigs(outDirCoverage, designDF, coverageDict[mitoName])

    ### process heteroplasmy tsv file
    # if there aren't heteroplasmy files, exit
    inVarTSV = inputDir / "out.heteroplasmy.tsv.gz"
    if(not inVarTSV.is_file()):
        exit()
    # if allowing skpping regen and the file already exists then read it in
    outDirVariants = Path(outputDir / "variants")
    if not outDirVariants.is_dir():
        outDirVariants.mkdir()
    
    outVarIntermediatePath = outDirVariants / "heteroplasmy.intermediate.tsv"
    if ((not outVarIntermediatePath.is_file()) or (args.regenfiles)):

        raw_variants = parseVariantTSV(inVarTSV)
        raw_variants['Genotype'] = raw_variants['Genotype'].replace({"0/0": "Ref",
                                "0/1": "Heteroplasmy",
                                "1/0": "Heteroplasmy", 
                                "1/1": "Homoplasmy",
                                "./.": "Missing"})
        annotationDF = getAnnotationsAll(raw_variants.query('(ID != "NORSID") and (Genotype != "Ref")')['ID'])
        merged_variants = raw_variants.reset_index().merge(designDF, how="left", on="Sample")
        merged_variants = merged_variants.merge(annotationDF, how="left", on="ID")
        merged_variants.to_csv(outVarIntermediatePath, sep="\t")
    else:
        merged_variants = pd.read_csv(outVarIntermediatePath, sep="\t")
    # get variant counts
    variants_only = merged_variants.query('Genotype != "Ref"')
    varCounts = pd.pivot_table(raw_variants, index="Sample", columns="Genotype", values="ADVar", aggfunc="count")
    varCounts = varCounts.merge(designDF, how="left", left_index=True, right_on="Sample")
    # write to heteroplasmy workbook
    
    # defined Excel workbook name earlier
    writer = pd.ExcelWriter(outDirWorkbooks / outExcelNameHeteroplasmy, engine='xlsxwriter')
    varCounts.to_excel(writer, "Summary Counts", index=False)
    for name,group in variants_only.groupby("Sample"):
        print(name)
        group.to_excel(writer, name, index=False)
    writer.save()
