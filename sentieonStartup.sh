#!/bin/bash
set -e

#run docker container

echo "Authenticating Sentieon License"
SENTIEON_KEY=JTTIF-HRBPG-GSSAC-RDGMS-KZMHL
credentials_file=~/credentials.json
project_file=${credentials_file}.project
python opt/gen_credentials.py $credentials_file "$SENTIEON_KEY" &

sleep 5
export SENTIEON_AUTH_MECH=proxy_GOOGLE
export SENTIEON_AUTH_DATA="$credentials_file"
read -r SENTIEON_JOB_TAG < "$project_file"
export SENTIEON_JOB_TAG

## Test DNS resolution before falling back to hard-coded IP
if nslookup gcp.sentieon.com > /dev/null; then
   lic_srvr_addr=gcp.sentieon.com
else
   lic_srvr_addr=35.188.47.218
 fi

export SENTIEON_LICENSE=$lic_srvr_addr:9003

#exec $0 "$@"

echo "Done Authenticating"
