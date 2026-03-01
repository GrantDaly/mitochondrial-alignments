# mitochondrial-alignments

This is a repository with wdl workflows and Python programs to align cell-free mitochondrial DNA in the fastq format.

The software then compiles statistics regarding coverage and insert size of the mitochondrial DNA and nuclear "NUMT" DNA.


# clone environment
conda env create --prefix envs -f conda-env-mito.yml

#build the docker image
docker build -t gdaly9000/mitochondrial .
docker push gdaly9000/mitochondrial
# verified versions of relevant software (Dockerhub sha256:851496157946b59912c7c80368f915c16e9011cf94914b8b84d649a16956178b from Jan 07, 2024 - Dec 22, 2025).
bwa                       0.7.17               h5bf99c6_8    bioconda
cutadapt                  4.6              py39hf95cd2a_1    bioconda
samtools 1.19-3-g62195d3
bcftools 1.12-34-gb1d9337
Using htslib 1.19-3-g67f3ab0f


# run tests for workbook generation
 python generate-workbooks.py \
-p tests/parameters.json \
-m tests/design.tsv \
-i tests/raw-postprocess/ -o tests/post-test/ \
-n test --regenfiles
