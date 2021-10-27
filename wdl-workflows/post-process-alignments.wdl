version 1.0

struct Alignment {
   String name
   File alignment
   File index
}

struct Bed {
   File bedFile
   String name
}

workflow coverageAndInsert {

   input {
    Array[File] sample_names
    Array[File] bam_files
    Array[File] bam_indexes
    
    Array[File] bed_names
    Array[File] bed_files
    
    File fasta
    File fasta_index
    }
    output {
    Array[Array[Array[File]]] output_files = process_beds.output_tsvs
    #Array[Array[File]] output_tsvs = combine.output_tsvs
    

    }
    Array[Int] sample_indices = range(length(sample_names))
    scatter (sample_index in sample_indices) {
    Array[Int] bed_indices = range(length(bed_names))    
    scatter (bed_index in bed_indices) {

        call Postprocess as process_beds{
           input:
               sample_name = sample_names[sample_index],
               sample_bam = bam_files[sample_index],
               sample_index  = bam_indexes[sample_index],
               groupName = bed_names[bed_index],
               bed = bed_files[bed_index]
        }
    }
    }
    
    
    }


task Postprocess {
  input {
    String sample_name
    File sample_bam
    File sample_index
    String groupName
    File bed
  }
  output {
    Array[File] output_tsvs = glob("*/*.tsv")
    #File cov_stats = "~{groupName}.coverage.stats.tsv"
    #File coverage = "~{groupName}.coverage.tsv"
    #File insert_stats = "~{groupName}.insert.stats.tsv"
    #File insert_hist = "~{groupName}.insert.hist.tsv"
  }
   runtime {
   docker: "gdaly9000/mitochondrial"
   cpu: 2
   memory: "16 GB"
   disks: "local-disk 200 SSD"
 }
  
  command {
    set -e
    echo "Post process Checksum " $(md5sum $(which inserts-and-cov))

    inserts-and-cov --bed ~{bed} --prefix ~{groupName} --bam ~{sample_bam} --outDir $(pwd) --name ~{sample_name}
  }
}

task RawVCF {
   input {
      String sample_name
      File sample_bam
      File sample_index
      File fasta
      File fasta_index
   }
   output {
      File vcf = "~{sample_name}.raw.vcf.gz"
      File index = "~{sample_name}.raw.vcf.gz.csi"
   }
   runtime {
   docker: "gdaly9000/mitochondrial"
   cpu: 4
   memory: "32 GB"
   disks: "local-disk 200 SSD"
 }
 command {
   bcftools --version
   bcftools mpileup -f ~{fasta} -r chrM -d 10000 -q 20 -Q 20 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR ~{sample_bam} | bcftools call -m -A | bcftools filter -e 'ALT=="."' -O z -o "~{sample_name}.raw.vcf.gz"
   bcftools index "~{sample_name}.raw.vcf.gz"
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

 
