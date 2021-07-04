To make the input files for the post processing
gsutil ls gs://south-al-genomics/novogene-trauma-4-20-21/alignments | awk '{print "{\"sample\":\""$1"\"},"}'
copy/paste into 
{"sampleList":[]}
