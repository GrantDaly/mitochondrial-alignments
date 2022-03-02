# mitochondrial-alignments

This is a repository with wdl workflows and Python programs to align cell-free mitochondrial DNA in the fastq format.

The software then compiles statistics regarding coverage and insert size of the mitochondrial DNA and nuclear "NUMT" DNA.


# clone environment
conda env create --prefix envs -f conda-env-mito.yml

#build the docker image
docker build -t gdaly9000/mitochondrial .
docker push gdaly9000/mitochondrial

# run tests for workbook generation
 python generate-workbooks.py \
-p tests/post-test/parameters.json \
-m tests/post-test/design.tsv \
-i tests/raw-postprocess/ -o tests/post-test/ \
-n test --regenfiles
