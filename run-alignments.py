import sys
from googleapiclient.discovery import build
import requests
import json
from pprint import pprint
import time
import glob
import pdb
import functools
import itertools


import re
from google.cloud import storage

from google.oauth2 import service_account

#credentials
credentials = service_account.Credentials.from_service_account_file(
    '/home/gtd1521_jagmail_southalabama_edu/huvecs/development/google-cloud/alignments-65005-fdd2a449ed54.json')

#could use this for web app
#credentials = GoogleCredentials.get_application_default()


scoped_credentials = credentials.with_scopes(
    ['https://www.googleapis.com/auth/cloud-platform'])
#build the service from api
service = build('lifesciences', 'v2beta', credentials=scoped_credentials)

#storage client
storageClient = storage.Client(project="alignments-65005", credentials=credentials)



# dryRun will avoid making alignments to test paramters
dryRun = False
inFileName = sys.argv[1]
print("loading input parameter json")
with open(inFileName) as f:
    sampleDict = json.load(f)

print("opening alignment pipeline json")
with open("alignment.pipeline") as alignPipeline:
    align_request_body = json.load(alignPipeline)
print("opening merge pipeline json")
with open("merge.pipeline") as mergePipeline:
    merge_bams_body = json.load(mergePipeline)

bucket = sampleDict['bucket']
bucket_object = storageClient.bucket(bucket)
parent = 'projects/alignments-65005/locations/us-central1'
prefixDirBlob =  sampleDict['blob-prefix']
prefixDir = "gs://" + bucket + "/" + prefixDirBlob

environment = align_request_body['pipeline']['environment']
environment['OUTPREFIX'] = prefixDir
environment['REFERENCEDIR'] = sampleDict['reference-dir']
environment['REFERENCE'] = sampleDict['reference']


sampleList = sampleDict['sampleList']
for sample in sampleList:
    

    sampleName = sample['sample']
    libraryName = sample['library']
    print("Sample %s" % sampleName)

    # check for existing files
    skipSample = False
    #if there are already sorted bams, delete them
    blobPrefix = prefixDirBlob + "/" + sampleName
    blobs = storageClient.list_blobs(bucket, prefix=blobPrefix)
    prog = re.compile(r'.sorted.bam')
    for blob in blobs:
        result = prog.search(blob.name)
        if(result):
            print("deleting %s" % blob.name)
            blob.delete()
    blobs = storageClient.list_blobs(bucket, prefix=blobPrefix)
    prog = re.compile(r'.sorted.merged.bam')
    #see if the merged entry is already there. If so, skip this sample
    for blob in blobs:
        result = prog.search(blob.name)
        if(result):
            print("existing merged bam %s , skipping to next sample" % blob.name)
            skipSample = True

    #don't want to repeat work so skipping if we already have the merged entry
    if(skipSample == True):
        continue

    
    numberLanes = len(sample['flowcellList'])
    for flowcell in sample['flowcellList']:
        flowcellName = flowcell['flowcell']
        laneName = flowcell['lane']
        fastq1 = flowcell['fastq1']
        fastq2 = flowcell['fastq2']

        print("Aligning Flowcell %s" % flowcellName)
        pipeline = align_request_body['pipeline']['environment']
        pipeline['FASTQ1'] = fastq1
        pipeline['FASTQ2'] = fastq2
        pipeline['SAMPLE'] = sampleName
        pipeline['FLOWCELL'] = flowcellName

        readgroup = r"@RG\tID:" + sampleName + "." + flowcellName + "." + laneName + r"\tSM:" + sampleName + r"\tLB:" + libraryName + r"\tPL:ILLUMINA"
        pipeline['READGROUP'] = readgroup

        pprint(pipeline)

        if not dryRun:
            workRequest = service.projects().locations().pipelines().run(parent=parent, body=align_request_body)
            workResponse = workRequest.execute()

            statusParent = workResponse['name']

            print("status url %s" % statusParent)
            notDone = True
            while notDone:
                statusRequest = service.projects().locations().operations().get(name=statusParent)
                statusResponse = statusRequest.execute()

                # if("events" in statusResponse['metadata'].keys()):
                #     machineList = statusResponse['metadata']['events']
                #     pprint(machineList)

                if("done" in statusResponse.keys()):
                    print("Done with " + sampleName + " " + flowcellName)
                    notDone = False
                else:
                    print("sleep 10 minutes")
                    time.sleep(600)
            
    if(numberLanes > 1):
        #merge alignments and delete unmerged
        mergeEnvs = merge_bams_body['pipeline']['environment']
        mergeEnvs['BAMDIR'] = prefixDir
        mergeEnvs['SAMPLE'] = sampleName
        #pdb.set_trace()
        if not dryRun:
            workRequest = service.projects().locations().pipelines().run(parent=parent, body=merge_bams_body)
            workResponse = workRequest.execute()

            statusParent = workResponse['name']
            notDone = True
            while notDone:
                statusRequest = service.projects().locations().operations().get(name=statusParent)
                statusResponse = statusRequest.execute()

                if("done" in statusResponse.keys()):
                    print("Done")
                    notDone = False
                else:
                    print("sleep 5 minutes")
                    time.sleep(300)
        
        print("Delete non-merged bams")
        blobPrefix = prefixDirBlob + "/" + sampleName
        blobs = storageClient.list_blobs(bucket, prefix=blobPrefix)
        prog = re.compile(r'.sorted.bam')
        for blob in blobs:
            print("deleting %s" % blob.name)
            result = prog.search(blob.name)
            if(result and not dryRun):
                blob.delete()
        print("Done Sample " + sampleName)

    elif(numberLanes == 1):
        
        blobPrefix = prefixDirBlob + "/" + sampleName
        newBamName = blobPrefix + "/" + sampleName + ".sorted.merged.bam"
        newBaiName = newBamName + ".bai"




        # copy bam
        blobs = storageClient.list_blobs(bucket, prefix=blobPrefix)

        prog = re.compile(r'.sorted.bam$')
        # copy bam
        for blob in blobs:
            #print("changing names from %s" % blob.name)
            result = prog.findall(blob.name)
            if(result and not dryRun):
                #blob.delete()
                # .sorted.merged.bam
                print("could copyt %s" % blob.name)
                print("new name  %s" % newBamName)
                _ = bucket_object.copy_blob(blob, bucket_object, new_name=newBamName)
                blob.delete()
        print("Done Sample " + sampleName)

        # copy bai
        blobs = storageClient.list_blobs(bucket, prefix=blobPrefix)
        prog = re.compile(r'.sorted.bam.bai')
        # copy bam
        for blob in blobs:
            #print("changing names from %s" % blob.name)
            result = prog.findall(blob.name)
            if(result and not dryRun):
                #blob.delete()
                # .sorted.merged.bam
                print("could copyt %s" % blob.name)
                print("new name  %s" % newBamName)
                _ = bucket_object.copy_blob(blob, bucket_object, new_name=newBaiName)
                blob.delete()
        print("Done Sample " + sampleName)

