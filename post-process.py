import sys
import json
import re
import pandas as pd
import numpy as np
from pathlib import Path
import pdb
import subprocess
import pyranges as pr

import pysam
from google.cloud import storage
from google.oauth2 import service_account

def getDepthSample(inBed, bamName, sample):
    inBed = pd.read_csv(inBed, sep="\t", names=["chr", "start", "end", "name", "score", "strand"])
    outDFList = []
    for bedCount, row in inBed.iterrows():
        chrom = row['chr']
        start = row['start']
        end = row['end']
        name = row['name']
        score = row['score']
        strand = row['strand']
        
        depthF1R2, depthF2R1 = getDepth(chrom, start, end, bamName)
        outDepthDFList = []
        for basePos, bases in enumerate(zip(depthF1R2, depthF2R1)):
            #print(basePos)
            countF1R2, countF2R1 = bases
            totalDepth = countF1R2 + countF2R1
            outDepthDFList.append({"Chromosome": chrom, "Start": start, "End": end,
                       "Name": name, "Score": score, "Strand": strand,
                        "Sample": sample, "Offset": basePos,
                                   "Forward_Depth_Mito": countF1R2, "Reverse_Depth_Mito": countF2R1,
                                   "Depth": totalDepth})
        outDFList.append(pd.DataFrame(outDepthDFList))
    outDF = pd.concat(outDFList, axis=0)
    return outDF

def readFilter(read):
    pair = read.is_proper_pair
    qual = read.mapping_quality

    if( pair and (qual >= 20)):
        return True
    else:
        return False
    
def readFilterF1R2(read):
    pair = read.is_proper_pair
    qual = read.mapping_quality

    
    if(qual >= 20):
        
        if((read.is_paired and read.is_proper_pair and read.mate_is_reverse and read.is_read1) or
           (read.is_paired and read.is_proper_pair and read.is_reverse and read.is_read2)):
            return True
    else:
        return False

def readFilterF2R1(read):
    pair = read.is_proper_pair
    qual = read.mapping_quality

    if(qual >= 20):
        if((read.is_paired and read.is_proper_pair and read.mate_is_reverse and read.is_read2) or
           (read.is_paired and read.is_proper_pair and read.is_reverse and read.is_read1)):
            return True
        
    else:
        return False

def getDepth(chrom, start, end, bamName, binSize=100):
    samfile = pysam.AlignmentFile(bamName, "rb")
    #print("before coverage pysam")
    
    coverageListF1R2 = []
    coverageListF2R1 = []
    
    tempCoverageListF1R2 = samfile.count_coverage(chrom, start=start, stop=end, read_callback=readFilterF1R2)
        
    tempCoverageListF1R2 = [sum(base) for base in zip(*tempCoverageListF1R2)]
    coverageListF1R2.extend(tempCoverageListF1R2)

    tempCoverageListF2R1 = samfile.count_coverage(chrom, start=start, stop=end, read_callback=readFilterF2R1)
        
    tempCoverageListF2R1 = [sum(base) for base in zip(*tempCoverageListF2R1)]
    coverageListF2R1.extend(tempCoverageListF2R1)
        
    samfile.close()

    return coverageListF1R2, coverageListF2R1

def getInsertsSample(inBed, bamName, sample, maxInsert):
    inBed = pd.read_csv(inBed, sep="\t", names=["chr", "start", "end", "name", "score", "strand"])
    outDFList = []
    for bedCount, row in inBed.iterrows():
        chrom = row['chr']
        start = row['start']
        end = row['end']        
        insertDF = pd.DataFrame(getInsertsByRegion(chrom, start, end, bamName, maxInsert))
        #         outInsertDFList = []
#         for index, insertCount in insertArray.iteritems():
#             outInsertDFList.append({"Sample": sample})
        outDFList.append(insertDF)
    outDF = pd.concat(outDFList, axis=0)
    return outDF

        
def getInsertsByRegion(chrom, start, end, bamName, maxInsert):
    rawList = []
    samfile = pysam.AlignmentFile(bamName, "rb")
    for read in samfile.fetch(chrom, start, end):
        #pair = read.is_proper_pair
        qual = read.mapping_quality

        insertSize = abs(read.template_length)

        # make sure there isn't excessive soft clipping. Usually occurs at insertion site
        # seems to return erroneous 0 length reads.
        # now an issue where there are genuinely clipped reads. also the idea that spikes were all from insertion events didn't pan out. smoothing function seemed to fix this
        readLength = read.infer_query_length()
        if not readLength:
            continue
        cigarStats, cigarBlocks  = read.get_cigar_stats()
        
        # numberMatches = cigarStats[0]
        # try:
        #     fracMatches = numberMatches / readLength
        # except TypeError:
        #      continue
        
        # if(fracMatches < 0.90):
        #     continue
        
        if(insertSize <= maxInsert):
            rawList.append(insertSize)
        
    samfile.close()
    return pd.Series(rawList, name="Fragments",dtype="int64")

def smoothDensity(inDF, binSize=3):
    inDFSorted = inDF.sort_values("Intervals")
    numEntries = len(inDFSorted)
    densityArray = inDF['% Density'].to_numpy()
    outArray = np.zeros(len(densityArray))
    centerOffset = binSize // 2

    for i in range(numEntries - binSize + 1):
#         print(densityArray[i: i + 3])
        tempMedian = np.median(densityArray[i: i + binSize])
        outArray[i + centerOffset] = tempMedian
    outSeries = pd.Series(outArray)
    return outSeries


def driveBedInsert(inBed, inBam, sample):

    maxInsert = 1000
    tempSample = getInsertsSample(inBed, inBam, sample, maxInsert)
    sampleStats = tempSample['Fragments'].describe()
    sampleStats['Sample'] = sample
        
    #histograms
    intervalLength = 1000
    sampleDens, sampleIntervals = np.histogram(tempSample['Fragments'],bins=np.arange(1,intervalLength + 2,1), density=True)
    sampleDens = pd.Series(sampleDens, name="Raw Density",dtype="float64")
    sampleIntervals = pd.Series(sampleIntervals[:intervalLength], name="Intervals",dtype="int64")
    sampleHist = pd.concat([sampleDens, sampleIntervals], axis=1)

    sampleHist['Sample'] = sample
    sampleHist['% Density'] = sampleHist['Raw Density'] * 100
    sampleHist['Smoothed % Density'] = smoothDensity(sampleHist,binSize=5)
    
    return sampleStats, sampleHist

credentials = service_account.Credentials.from_service_account_file(
    '/home/gtd1521_jagmail_southalabama_edu/huvecs/development/google-cloud/alignments-65005-fdd2a449ed54.json')

#could use this for web app
#credentials = GoogleCredentials.get_application_default()

scoped_credentials = credentials.with_scopes(
    ['https://www.googleapis.com/auth/cloud-platform'])
#build the service from api
#service = build('lifesciences', 'v2beta', credentials=scoped_credentials)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--name', default="DAMPs")
parser.add_argument('inFile')
parser.add_argument('workingDir')
parser.add_argument("--noRegen", help="if the sample was already processed, can skip",
                    action="store_true")
parser.add_argument("--datetime",action="store_true")
args = parser.parse_args()
# from bottom need to update
outDir = Path(args.workingDir)
if not outDir.is_dir():
    print("Need existing directory with design.csv file")
    exit()


workingDir = outDir
outName = args.name
inFileName = args.inFile



dryRun = False

with open(inFileName) as f:
    sampleDict = json.load(f)



#storage client
storageClient = storage.Client(project="alignments-65005", credentials=credentials)

#bucket = sampleDict['bucket']
#parent = 'projects/alignments-65005/locations/us-central1'
#prefixDirBlob =  sampleDict['blob-prefix']
#prefixDir = "gs://" + bucket + "/" + prefixDirBlob

#outDirVariants = Path("variants")
#if not outDirVariants.is_dir():
#    outDirVariants.mkdir()

# Depth lists for each bed file
mitoDepthList = []
numtDepthList = []
dayamaDepthList = []
shuffleDepthList = []

# Inserts lists for each bed file
mitoInsertStatsList = []
numtInsertStatsList = []
dayamaInsertStatsList = []
shuffleInsertStatsList = []

mitoInsertHistList = []
numtInsertHistList = []
dayamaInsertHistList = []
shuffleInsertHistList = []

# output directory prefix is second command-line argument

outDirCov = outDir / "coverage-outputs"
if not outDirCov.is_dir():
    outDirCov.mkdir()

outDirInserts = outDir / "insert-outputs"
if not outDirInserts.is_dir():
    outDirInserts.mkdir()

tempDir = outDir / "temp"
if not tempDir.is_dir():
    tempDir.mkdir()


sampleList = sampleDict['sampleList']
for sample in sampleList:
    sampleURL = sample['sample']

    # in json files need trailing "/"
    sampleName = sampleURL.split("/")[-2]
    
    print(sampleName)

    # if in regen mode, check if files exist and skip ifso    
    numtDepthCovName = outDirCov / (sampleName + ".numt.coverage.tsv")
    tempInsertStatsNumtName = outDirInserts / (sampleName + ".numt.insert.stats.tsv")
    tempInsertHistNumtName = outDirInserts / (sampleName + ".numt.insert.hist.tsv")
    shufDepthCovName = outDirCov / (sampleName + ".shuf.coverage.tsv")
    dayamaDepthCovName = outDirCov / (sampleName + ".dayama.coverage.tsv")
    mitoDepthCovName = outDirCov / (sampleName + ".mito.coverage.tsv")
    insertHistMitoName = outDirInserts / (sampleName + ".mito.insert.hist.tsv")
    tempInsertStatsMitoName = outDirInserts / (sampleName + ".mito.insert.stats.tsv")

    
    if(args.noRegen):
        
        if(numtDepthCovName.exists() and
                tempInsertStatsNumtName.exists() and
                tempInsertStatsNumtName.exists() and
                shufDepthCovName.exists() and 
                dayamaDepthCovName.exists() and
                mitoDepthCovName.exists() and
                insertHistMitoName.exists() and
                tempInsertStatsMitoName.exists()):
            continue
    

    sampleRegex = re.compile(r'gs:\/\/([a-zA-Z0-9_-]+)(\/)(.+)')
    sampleMatch = sampleRegex.match(sampleURL)

    bucket = sampleMatch.group(1)
    blobPrefix = sampleMatch.group(3)
    
    progBam = re.compile(r'.sorted.merged.bam$')
    progBai = re.compile(r'.sorted.merged.bam.bai$')

    # if the data were sequenced on only one lane they won't be merged
    
    backupBam = re.compile(r'.bam$')
    backupBai = re.compile(r'.bam.bai$')


    tempBamName = str( tempDir / "temp.bam")
    tempBaiName = tempBamName + ".bai"

    hadBam = False
    blobs = storageClient.list_blobs(bucket, prefix=blobPrefix)
    for blob in blobs:
        
        result = progBam.search(blob.name)

        if not result:
            result = backupBam.search(blob.name)
        if(result):
            hadBam = True
            print("Downloading %s" % blob.name)
            blob.download_to_filename(tempBamName)
    if not hadBam:
        print("No bam")
        quit()

    hadBai = False
    blobs = storageClient.list_blobs(bucket, prefix=blobPrefix)
    for blob in blobs:
        
        result = progBai.search(blob.name)

        if not result:
            result = backupBai.search(blob.name)
        if(result):
            hadBai = True
            print("Downloading %s" % blob.name)
            blob.download_to_filename(tempBaiName)

    if not hadBai:
        print("No bam index")
        quit()

        #tempMito = getDepthSample(inBed, inBam, sample)
        #mitoList.append(tempMito)

        print("Generating Coverages")

    
    print("Calculating Numt Stats")
    numtBedName = "beds/numt-GRCh38.fixed.bed"
    tempNumtDepth = getDepthSample(numtBedName,tempBamName , sampleName)

    #numtDepthList.append(tempNumtDepth)
    
    tempNumtDepth.to_csv(numtDepthCovName, sep="\t", index=False)
    
    tempInsertStatsNumt, tempInsertHistNumt = driveBedInsert(numtBedName, tempBamName, sampleName)

    # this is small enough that we can keep in memory
    #numtInsertStatsList.append(tempInsertStatsNumt)
    #numtInsertHistList.append(tempInsertHistNumt)

    
    pd.DataFrame(tempInsertStatsNumt).T.to_csv(tempInsertStatsNumtName, sep="\t")

    tempInsertHistNumt.to_csv(tempInsertHistNumtName, sep="\t",index=None)


    print("Calculating Shuffled Stats")
    shufBedName = "beds/shuf-numt.bed"
    tempShufDepth = getDepthSample(shufBedName,tempBamName , sampleName)

    tempShufDepth.to_csv(shufDepthCovName, sep="\t", index=False)
    #shuffleDepthList.append(tempShufDepth)

    # don't want to keep shuffled inserts because sparse
    # tempInsertStatsShuf, tempInsertHistShuf = driveBedInsert(shufBedName, tempBamName, sampleName)
    # shuffleInsertStatsList.append(tempInsertStatsShuf)
    # shuffleInsertHistList.append(tempInsertHistShuf)
    
    print("Calculating Dayama Stats")
    dayamaBedName = "beds/dayama.grch38.200slop.bed"
    tempDayamaDepth = getDepthSample(dayamaBedName,tempBamName , sampleName)
    #dayamaDepthList.append(tempDayamaDepth)
    tempDayamaDepth.to_csv(dayamaDepthCovName, sep="\t", index=False)


    # also sparse so don't want inserts
    # tempInsertStatsDayama, tempInsertHistDayama = driveBedInsert(dayamaBedName, tempBamName, sampleName)
    # dayamaInsertStatsList.append(tempInsertStatsDayama)
    # dayamaInsertHistList.append(tempInsertHistDayama)


    print("Calculating Mito Stats")
    mitoBedName = "beds/mito.full.human.bed"
    print("Coverage")
    
    tempMitoDepth = getDepthSample(mitoBedName,tempBamName , sampleName)

    tempMitoDepth.to_csv(mitoDepthCovName , sep="\t", index=False)
    #mitoDepthList.append(tempMitoDepth)
    print("Done with Coverage")

    print("Inserts")
    tempInsertStatsMito, tempInsertHistMito = driveBedInsert(mitoBedName, tempBamName, sampleName)
    print("Done with Inserts")
    # this is small enough that we can keep in memory
    #mitoInsertStatsList.append(tempInsertStatsMito)
    #numtInsertHistList.append(tempInsertHistNumt)

    tempInsertHistMito.to_csv(insertHistMitoName , sep="\t", index=False)
    
    pd.DataFrame(tempInsertStatsMito).T.to_csv(tempInsertStatsMitoName, sep="\t")
    
    

    # prep variants
    # print("running variant calls")
    # subprocess.run(["./run-variant-calls.sh", sampleName, tempBamName])
    Path(tempBamName).unlink()
    Path(tempBaiName).unlink()

    print("Done with Sample %s" % sampleName)

print("Done with processing")
# concat Depth Samples

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



# Depths
# aggregating the numt and shuffled numt because they are the largest. should help avoid out of memory exception
numtDepthList = []
for depthFile in outDirCov.glob("*.numt.coverage.tsv"):
    tempCov = pd.read_csv(depthFile, sep="\t",dtype=covDtypes)
    tempCov = tempCov.groupby(["Name", "Chromosome", "Start", "End", "Strand","Score","Sample"])["Depth"].agg('sum').reset_index()
    numtDepthList.append(tempCov)

concatDepthNumt = pd.concat(numtDepthList, axis=0)
concatDepthNumt.to_csv(outDirCov / "Numt-Coverages-Agg.tsv", sep="\t", index=False)
del numtDepthList

# Depth lists for each bed file
shuffleDepthList = []

for depthFile in outDirCov.glob("*.shuf.coverage.tsv"):
    tempCov = pd.read_csv(depthFile, sep="\t",dtype=covDtypes)
    tempCov = tempCov.groupby(["Name", "Chromosome", "Start", "End", "Strand","Score","Sample"])["Depth"].agg('sum').reset_index()
    shuffleDepthList.append(tempCov)
concatDepthShuf = pd.concat(shuffleDepthList, axis=0)
concatDepthShuf.to_csv(outDirCov / "Shuffled-Coverages-Agg.tsv", sep="\t", index=False)
del concatDepthShuf

dayamaDepthList = []
for depthFile in outDirCov.glob("*.shuf.coverage.tsv"):
    tempCov = pd.read_csv(depthFile, sep="\t",dtype=covDtypes)
    dayamaDepthList.append(tempCov)
concatDepthDayama = pd.concat(dayamaDepthList, axis=0)
concatDepthDayama.to_csv(outDirCov / "Dayama-Coverages.tsv", sep="\t", index=False)

mitoDepthList = []
for depthFile in outDirCov.glob("*.mito.coverage.tsv"):
    tempCov = pd.read_csv(depthFile, sep="\t",dtype=covDtypes)
    mitoDepthList.append(tempCov)

concatDepthMito = pd.concat(mitoDepthList, axis=0)
concatDepthMito.to_csv(outDirCov / "Mito-Coverages.tsv", sep="\t", index=False)


# concat Inserts

# Inserts lists for each bed file

print("Writing intermediate tsv files")
numtInsertStatsList = []
for insertFile in outDirInserts.glob("*.numt.insert.stats.tsv"):
    tempInsert = pd.read_csv(insertFile, sep="\t")
    numtInsertStatsList.append(tempInsert)
concatInsertStatsNumt = pd.concat(numtInsertStatsList, axis=0) 
concatInsertStatsNumt.to_csv(outDirInserts / "Numt-InsertStats.tsv", sep="\t", index=False)

mitoInsertStatsList = []
for insertFile in outDirInserts.glob("*.mito.insert.stats.tsv"):
    tempInsert = pd.read_csv(insertFile, sep="\t")
    mitoInsertStatsList.append(tempInsert)
concatInsertStatsMito = pd.concat(mitoInsertStatsList, axis=0) 
concatInsertStatsMito.to_csv(outDirInserts / "Mito-InsertStats.tsv", sep="\t", index=False)


numtInsertHistList = []
for insertFile in outDirInserts.glob("*.numt.insert.hist.tsv"):
    tempInsert = pd.read_csv(insertFile, sep="\t")
    numtInsertHistList.append(tempInsert)
concatInsertHistNumt = pd.concat(numtInsertHistList, axis=0)
concatInsertHistNumt.to_csv(outDirInserts / "Numt-InsertHist.tsv", sep="\t", index=False)

mitoInsertHistList = []
for insertFile in outDirInserts.glob("*.mito.insert.hist.tsv"):
    tempInsert = pd.read_csv(insertFile, sep="\t")
    mitoInsertHistList.append(tempInsert)
concatInsertHistMito = pd.concat(mitoInsertHistList, axis=0)
concatInsertHistMito.to_csv(outDirInserts / "Mito-InsertHist.tsv", sep="\t", index=False)








# these were the ones kept in memory
# concatInsertStatsNumt = pd.DataFrame(numtInsertStatsList).reset_index(drop=True)
# #concatInsertStatsShuf = pd.DataFrame(shuffleInsertStatsList).reset_index(drop=True)
# #concatInsertStatsDayama = pd.DataFrame(dayamaInsertStatsList).reset_index(drop=True)
# concatInsertStatsMito = pd.DataFrame(mitoInsertStatsList).reset_index(drop=True)

# concatInsertHistNumt = pd.concat(numtInsertHistList, axis=0)
# #concatInsertHistShuf = pd.concat(shuffleInsertHistList, axis=0)
# #concatInsertHistDayama = pd.concat(dayamaInsertHistList, axis=0)
# concatInsertHistMito = pd.concat(mitoInsertHistList, axis=0)



# workbook stage

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

numtCovCutdown = pd.read_csv(coverageDir / "Numt-Coverages-Agg.tsv", sep="\t",dtype=covDtypes)

shufNumtCovCutdown = pd.read_csv(coverageDir / "Shuffled-Coverages-Agg.tsv", sep="\t",dtype=covDtypes)

dayamaCovCutdown = pd.read_csv(coverageDir / "Dayama-Coverages.tsv", sep="\t",dtype=covDtypes)
# these are the coverages aligned only to the mitochondria


#aggregate by numt
#covByNumt = numtCovCutdown.groupby(["Name", "Chromosome", "Start", "End", "Strand","Score","Sample"])["Depth"].agg('sum').reset_index()
covByNumt = pd.read_csv(coverageDir / "Numt-Coverages-Agg.tsv", sep="\t",dtype=covDtypes)
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

now = datetime.datetime.now()
if(args.datetime):
    writer = pd.ExcelWriter(outDirWorkbooks / (outName + str(now) + '.xlsx'), engine='xlsxwriter')
else:
    writer = pd.ExcelWriter(outDirWorkbooks / (outName + '.xlsx'), engine='xlsxwriter')

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

rawDesignList = []
normDesignList = []
for name, group in humanMtNorm.groupby("Sample"):

    designEntry = designDF.query('Sample == @name')
    patient = designEntry['Patient ID'].to_list()[0]
    sampleName = designEntry['Sample'].to_list()[0]
    condition = designEntry['Condition'].to_list()[0]
    
    outGroup = group.copy().reset_index(drop=True)
    outGroup['Start'] = outGroup['Offset']
    outGroup['End']  = outGroup['Start'] + 1
    outGroupNorm = outGroup.copy()

#    outGroup['Count'] = outGroup['Depth_Mito'] 
#    outGroupNorm['Count'] = outGroup['Norm Depth']

    outNameRaw = name + "-raw"
    outNameNorm = name + "-norm"

#     filenameList.append(outName)

    # only outputting bigwigs for now.
    
    outGroupRawAll = outGroup.loc[:, ["Chromosome", "Start", "End", "Depth_Mito", "Forward_Depth_Mito", "Reverse_Depth_Mito"]]
    
    #outGroupRawAll['Percent_Diff_Strand'] = ((outGroupRawAll['Forward_Depth_Mito'] - outGroupRawAll['Reverse_Depth_Mito']) /  outGroupRawAll['Reverse_Depth_Mito']) * 100
    #pdb.set_trace()
    outGroupRawAll['LogFC_Diff_Strand'] = np.log2((outGroupRawAll['Forward_Depth_Mito'] + 1) /  (outGroupRawAll['Reverse_Depth_Mito'] + 1))
    
    #outGroup.to_csv(outDirRaw / (outNameRaw + ".bg"), sep="\t", index=False, header=False)
    
    outGroupNorm = outGroupNorm.loc[:, ["Chromosome", "Start", "End", "Norm Depth"]]
    #outGroupNorm.to_csv(outDirNorm / (outNameNorm + ".bg"), sep="\t", index=False, header=False)

    # pyranges to_bigwig not working
    tempPyRangesRaw = pr.PyRanges(outGroupRawAll)
    
    outRawString = str(outDirRaw / (outNameRaw + ".bw"))
    outRawStringFor = str(outDirRaw / (outNameRaw + ".F1R2.bw"))
    outRawStringRev = str(outDirRaw / (outNameRaw + ".F2R1.bw"))
    outRawStringDiff = str(outDirRaw / (outNameRaw + ".StrandDiff.bw"))
    tempPyRangesRaw.to_bigwig(path = outRawString, value_col="Depth_Mito")
    #pdb.set_trace()
    rawDesignList.append({'TRACK_ID': outNameRaw + ".bw", 'PARTICIPANT_ID': patient,
                          'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand': "Both",
                          'Condition-Strand': 'NA'})
    tempPyRangesRaw.to_bigwig(path = outRawStringFor, value_col="Forward_Depth_Mito")
    rawDesignList.append({'TRACK_ID': outNameRaw + ".F1R2.bw", 'PARTICIPANT_ID': patient,
                          'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand':"Forward",
                          "Condition-Strand": condition + "-Forward"})
 
    tempPyRangesRaw.to_bigwig(path = outRawStringRev, value_col="Reverse_Depth_Mito")
    rawDesignList.append({'TRACK_ID': outNameRaw + ".F2R1.bw", 'PARTICIPANT_ID': patient,
                          'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand':"Reverse",
                          'Condition-Strand': condition + "-Reverse"})

    tempPyRangesRaw.to_bigwig(path = outRawStringDiff, value_col="LogFC_Diff_Strand")
    rawDesignList.append({'TRACK_ID': outNameRaw + ".StrandDiff.bw", 'PARTICIPANT_ID': patient,
                          'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand':"LogFC",
                          'Condition-Strand': "NA"})
 
    
    tempPyRangesNorm = pr.PyRanges(outGroupNorm)
    outNormString = str(outDirNorm / (outNameNorm + ".bw"))
    tempPyRangesNorm.to_bigwig(path = outNormString, value_col="Norm Depth")
    normDesignList.append({'TRACK_ID': outNameNorm + ".bw", 'PARTICIPANT_ID': patient,
                          'SAMPLE_ID': sampleName, 'Condition': condition, 'Strand':"Both"})
 

outRawDesignDF = pd.DataFrame(rawDesignList)
outRawDesignDF['DATA_Type'] = "Expression"

outRawDesignDF.loc[:,["TRACK_ID", "DATA_Type", "PARTICIPANT_ID", "SAMPLE_ID", "Condition","Strand","Condition-Strand"]]\
                         .to_csv(outDirRaw / "raw-attributes.txt", sep="\t", index=None)


outNormDesignDF = pd.DataFrame(normDesignList)
outNormDesignDF['DATA_Type'] = "Expression"

outNormDesignDF.loc[:,["TRACK_ID", "DATA_Type", "PARTICIPANT_ID", "SAMPLE_ID", "Condition","Strand"]]\
                         .to_csv(outDirNorm / "norm-attributes.txt", sep="\t", index=None)
