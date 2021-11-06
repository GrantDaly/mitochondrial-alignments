version 1.0

workflow CallHeteroplasmies {

   input {
    Array[File] sample_names
    Array[File] bam_files
    Array[File] bam_indexes
    
    File fasta
    File fasta_index
    }
    output {

    File vcf = raw_vcf.vcf
    File vcf_index = raw_vcf.index
    
    File called_vcf = raw_vcf.called_vcf
    File called_vcf_index = raw_vcf.called_vcf_index

    
    }
    #Array[Int] sample_indices = range(length(sample_names))

        call GenerateRawVCF as raw_vcf{
           input:
               sample_bams = bam_files,
               sample_indexes  = bam_indexes,
               fasta = fasta,
               fasta_index = fasta_index
        }
    }



task GenerateRawVCF {
   input {
      Array[File] sample_bams
      Array[File] sample_indexes
      File fasta
      File fasta_index
   }
   output {

      File vcf = "out.raw.vcf.gz"
      File index = "out.raw.vcf.gz.csi"
      
      File called_vcf = "out.heteroplasmy.vcf.gz"
      File called_vcf_index = "out.heteroplasmy.vcf.gz.csi"

   }
   runtime {
   docker: "gdaly9000/mitochondrial"
   cpu: 4
   memory: "32 GB"
   disks: "local-disk 200 SSD"
 }
 command {
   bcftools --version
   bcftools mpileup -f ~{fasta} -r chrM -d 10000 -q 20 -Q 20 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR ~{sep=" " sample_bams}| bcftools call -m -A | bcftools filter -e 'ALT=="."' -O z -o "out.raw.vcf.gz"
   bcftools index "out.raw.vcf.gz"
   

   bcftools plugin heteroplasmy --version

   bcftools plugin heteroplasmy -O z -o "out.heteroplasmy.vcf.gz" "out.raw.vcf.gz"
   bcftools index "out.heteroplasmy.vcf.gz" 

 }
}

task ProcessRawVCFs {
   input {
      Array[File] input_vcfs
      Array[File] input_indexes
      
      File fasta
      File fasta_index
   }
   output {
      
      File out_vcf = "out.raw.vcf.gz"
      File out_vcf_index = "out.raw.vcf.gz.csi"
      
      File out_heteroplasmy = "out.heteroplasmy.annotated.vcf.gz"
      File out_heteroplasmy_index = "out.heteroplasmy.annotated.vcf.gz.csi"
      
      File out_heteroplasmy_tsv = "out.heteroplasmy.tsv.gz"
   }
   runtime {
   docker: "gdaly9000/mitochondrial"
   cpu: 4
   memory: "16 GB"
 }
 command <<<
   set -e
   set -o
   # merge into one vcf
   bcftools merge -O z -o "out.raw.vcf.gz" ~{sep=" " input_vcfs}
   bcftools index "out.raw.vcf.gz"
   
   # call heteroplasmy
   bcftools plugin heteroplasmy -O z -o "out.heteroplasmy.vcf.gz" "out.raw.vcf.gz"
   bcftools index "out.heteroplasmy.vcf.gz"
   
   # annotate
   zcat "out.heteroplasmy.vcf.gz" | sed 's/chrM/chrMT/g' > "variant.vcf"
   wget  "https://api.ncbi.nlm.nih.gov/variation/v0/vcf/file/set_rsids?assembly=GCF_000001405.38" --post-file "variant.vcf" -O "dbsnp.noheader.vcf"

   # concatenate header
   cat <( bcftools view -h "out.heteroplasmy.vcf.gz" ) <(sed 's/chrMT/chrM/g' "dbsnp.noheader.vcf") | bcftools view -O z -o "out.heteroplasmy.annotated.vcf.gz"
   bcftools index "out.heteroplasmy.annotated.vcf.gz"

   # create tsv file for downstream analysis
   bcftools norm -m - --fasta-ref "~{fasta}" "out.heteroplasmy.annotated.vcf.gz" | bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\tSample:ADRef:ADVar:ADFor-Ref:ADFor-Var:ADRev-Ref:ADRev-Var:SBFisher:ArtifactScore:MinorFreq:Genotype[\t%SAMPLE\t%AD{0}\t%AD{1}\t%ADF{0}\t%ADF{1}\t%ADR{0}\t%ADR{1}\t%SP\t%ART\t%AF\t%GT]\n" - | gzip > "out.heteroplasmy.tsv.gz"

 >>>
}

 
