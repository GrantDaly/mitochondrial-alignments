import sys
import json
import re
import pandas as pd
import numpy as np
from pathlib import Path
import pdb
import subprocess

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
        
        tempDepth = getDepth(chrom, start, end, bamName)
        outDepthDFList = []
        for basePos, base in enumerate(tempDepth):
            #print(basePos)
            outDepthDFList.append({"Chromosome": chrom, "Start": start, "End": end,
                       "Name": name, "Score": score, "Strand": strand,
                       "Sample": sample, "Offset": basePos, "Depth": base})
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
        
def getDepth(chrom, start, end, bamName, binSize=100):
    samfile = pysam.AlignmentFile(bamName, "rb")
    #print("before coverage pysam")
    binRanges = np.arange(start, end + 1, binSize)
    coverageList = []
    
    for tempStart in binRanges:
        tempEnd = min(tempStart + binSize, end)
        #print("Bin Starting with %d" % tempStart)
        if not (tempStart < tempEnd):
            continue
        tempCoverageList = samfile.count_coverage(chrom, start=tempStart, stop=tempEnd, read_callback=readFilter)
        
        tempCoverageList = [sum(base) for base in zip(*tempCoverageList)]
        coverageList.extend(tempCoverageList)
        
    samfile.close()

    return coverageList

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
        readLength = read.infer_query_length()
        if not readLength:
            continue
        cigarStats, cigarBlocks  = read.get_cigar_stats()
        
        numberMatches = cigarStats[0]
        try:
            fracMatches = numberMatches / readLength
        except TypeError:
             continue
        
        if(fracMatches < 0.90):
            continue
        
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

dryRun = False
inFileName = sys.argv[1]
with open(inFileName) as f:
    sampleDict = json.load(f)

#storage client
storageClient = storage.Client(project="alignments-65005", credentials=credentials)

#bucket = sampleDict['bucket']
#parent = 'projects/alignments-65005/locations/us-central1'
#prefixDirBlob =  sampleDict['blob-prefix']
#prefixDir = "gs://" + bucket + "/" + prefixDirBlob

outDirVariants = Path("variants")
if not outDirVariants.is_dir():
    outDirVariants.mkdir()

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


sampleList = sampleDict['sampleList']
for sample in sampleList:
    sampleURL = sample['sample']

    # in json files need trailing "/"
    sampleName = sampleURL.split("/")[-2]
    
    print(sampleName)

    sampleRegex = re.compile(r'gs:\/\/([a-zA-Z0-9_-]+)(\/)(.+)')
    sampleMatch = sampleRegex.match(sampleURL)

    bucket = sampleMatch.group(1)
    blobPrefix = sampleMatch.group(3)
    
    progBam = re.compile(r'.sorted.merged.bam$')
    progBai = re.compile(r'.sorted.merged.bam.bai$')

    # if the data were sequenced on only one lane they won't be merged
    
    backupBam = re.compile(r'.bam$')
    backupBai = re.compile(r'.bam.bai$')


    hadBam = False
    blobs = storageClient.list_blobs(bucket, prefix=blobPrefix)
    for blob in blobs:
        
        result = progBam.search(blob.name)

        if not result:
            result = backupBam.search(blob.name)
        if(result):
            hadBam = True
            print("Downloading %s" % blob.name)
            blob.download_to_filename("/mnt/disks/scratch/temp.bam")
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
            blob.download_to_filename("/mnt/disks/scratch/temp.bam.bai")

    if not hadBai:
        print("No bam index")
        quit()

        #tempMito = getDepthSample(inBed, inBam, sample)
        #mitoList.append(tempMito)

        print("Generating Coverages")

    tempBamName = "/mnt/disks/scratch/temp.bam"



    print("Calculating Numt Stats")
    numtBedName = "beds/numt-GRCh38.fixed.bed"
    tempNumtDepth = getDepthSample(numtBedName,tempBamName , sampleName)
    numtDepthList.append(tempNumtDepth)
    
    tempInsertStatsNumt, tempInsertHistNumt = driveBedInsert(numtBedName, tempBamName, sampleName)
    numtInsertStatsList.append(tempInsertStatsNumt)
    numtInsertHistList.append(tempInsertHistNumt)

    print("Calculating Shuffled Stats")
    shufBedName = "beds/shuf-numt.bed"
    tempShufDepth = getDepthSample(shufBedName,tempBamName , sampleName)
    shuffleDepthList.append(tempShufDepth)

    # don't want to keep shuffled inserts because sparse
    # tempInsertStatsShuf, tempInsertHistShuf = driveBedInsert(shufBedName, tempBamName, sampleName)
    # shuffleInsertStatsList.append(tempInsertStatsShuf)
    # shuffleInsertHistList.append(tempInsertHistShuf)
    
    print("Calculating Dayama Stats")
    dayamaBedName = "beds/dayama.grch38.200slop.bed"
    tempDayamaDepth = getDepthSample(dayamaBedName,tempBamName , sampleName)
    dayamaDepthList.append(tempDayamaDepth)


    # also sparse so don't want inserts
    # tempInsertStatsDayama, tempInsertHistDayama = driveBedInsert(dayamaBedName, tempBamName, sampleName)
    # dayamaInsertStatsList.append(tempInsertStatsDayama)
    # dayamaInsertHistList.append(tempInsertHistDayama)


    print("Calculating Mito Stats")
    mitoBedName = "beds/mito.full.human.bed"
    print("Coverage")
    tempMitoDepth = getDepthSample(mitoBedName,tempBamName , sampleName)
    mitoDepthList.append(tempMitoDepth)
    print("Done with Coverage")

    print("Inserts")
    tempInsertStatsMito, tempInsertHistMito = driveBedInsert(mitoBedName, tempBamName, sampleName)
    print("Done with Inserts")
    mitoInsertStatsList.append(tempInsertStatsMito)
    mitoInsertHistList.append(tempInsertHistMito)

    # prep variants
    # print("running variant calls")
    # subprocess.run(["./run-variant-calls.sh", sampleName, tempBamName])
    print("Done with Sample %s" % sampleName)

print("Done with processing")
# concat Depth Samples
concatDepthNumt = pd.concat(numtDepthList, axis=0)
concatDepthShuf = pd.concat(shuffleDepthList, axis=0)
concatDepthDayama = pd.concat(dayamaDepthList, axis=0)
concatDepthMito = pd.concat(mitoDepthList, axis=0)


# output directory prefix is second command-line argument
outDir = Path(sys.argv[2])
outDirCov = outDir / "coverage-outputs"
if not outDirCov.is_dir():
    outDirCov.mkdir()

concatDepthNumt.to_csv(outDirCov / "Numt-Coverages.tsv", sep="\t", index=False)
concatDepthShuf.to_csv(outDirCov / "Shuffled-Coverages.tsv", sep="\t", index=False)
concatDepthDayama.to_csv(outDirCov / "Dayama-Coverages.tsv", sep="\t", index=False)
concatDepthMito.to_csv(outDirCov / "Mito-Coverages.tsv", sep="\t", index=False)

# concat Inserts

concatInsertStatsNumt = pd.DataFrame(numtInsertStatsList).reset_index(drop=True)
#concatInsertStatsShuf = pd.DataFrame(shuffleInsertStatsList).reset_index(drop=True)
#concatInsertStatsDayama = pd.DataFrame(dayamaInsertStatsList).reset_index(drop=True)
concatInsertStatsMito = pd.DataFrame(mitoInsertStatsList).reset_index(drop=True)

concatInsertHistNumt = pd.concat(numtInsertHistList, axis=0)
#concatInsertHistShuf = pd.concat(shuffleInsertHistList, axis=0)
#concatInsertHistDayama = pd.concat(dayamaInsertHistList, axis=0)
concatInsertHistMito = pd.concat(mitoInsertHistList, axis=0)

outDirInserts = outDir / "insert-outputs"
if not outDirInserts.is_dir():
    outDirInserts.mkdir()

print("Writing intermediate tsv files")
concatInsertStatsNumt.to_csv(outDirInserts / "Numt-InsertStats.tsv", sep="\t", index=False)
#concatInsertStatsShuf.to_csv(outDirInserts / "Shuffled-InsertStats.tsv", sep="\t", index=False)
#concatInsertStatsDayama.to_csv(outDirInserts / "Dayama-InsertStats.tsv", sep="\t", index=False)
concatInsertStatsMito.to_csv(outDirInserts / "Mito-InsertStats.tsv", sep="\t", index=False)

concatInsertHistNumt.to_csv(outDirInserts / "Numt-InsertHist.tsv", sep="\t", index=False)
#concatInsertHistShuf.to_csv(outDirInserts / "Shuffled-InsertHist.tsv", sep="\t", index=False)
#concatInsertHistDayama.to_csv(outDirInserts / "Dayama-InsertHist.tsv", sep="\t", index=False)
concatInsertHistMito.to_csv(outDirInserts / "Mito-InsertHist.tsv", sep="\t", index=False)

