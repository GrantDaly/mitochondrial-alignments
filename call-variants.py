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
with open("variants.pipeline") as variantPipeline:
    variant_request_body = json.load(variantPipeline)

environment = variant_request_body['pipeline']['environment']
environment['FASTAPREFIX'] = sampleDict['FASTAPREFIX']
parent = 'projects/alignments-65005/locations/us-central1'

sampleList = sampleDict['sampleList']
for sample in sampleList:

    
    environment['INBAM'] = sample['INBAM']
    environment['OUTDIR'] = sample['OUTDIR']
    environment['OUTNAME'] = sample['OUTNAME']
    
        
    print("Calling  %s" % environment['INBAM'])

    

    pprint(environment)

    
    workRequest = service.projects().locations().pipelines().run(parent=parent, body=variant_request_body)
    workResponse = workRequest.execute()

    statusParent = workResponse['name']

    print("status url %s" % statusParent)
    notDone = True
    while notDone:
        statusRequest = service.projects().locations().operations().get(name=statusParent)
        statusResponse = statusRequest.execute()

        if("events" in statusResponse['metadata'].keys()):
             machineList = statusResponse['metadata']['events']
             pprint(machineList)

        if("done" in statusResponse.keys()):
            print("Done with " + environment['INBAM'])
            notDone = False
        else:
            print("sleep 1 minutes")
            time.sleep(60)
            
