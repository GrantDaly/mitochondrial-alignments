version 1.0

workflow alignSample {

   input {

    File ref_fasta
    File ref_fai
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    
   String sample_name
   #Array[String] library_names
   #Array[String] flowcell_names
   #Array[String] flowcell_lanes
   #Array[File] fastq_ones
   #Array[File] fastq_twos
   File fastq_tsv
    }

    output {
    File out_bam = merge.out_bam
    File out_bam_index = merge.out_bam_index
    Array[File] out_artifacts = Artifacts.out_artifacts
    File out_bam_stats = merge.bam_stats
    Int raw_read_count = merge.raw_read_count
    }

    #Array[Int] indices = range(length(fastq_ones))
    Array[Array[String]] fastq_tsv_array = read_tsv(fastq_tsv)
    scatter (flowcell in fastq_tsv_array) {
        
        call AlignAndSort as align{
           input:
               ref_fasta = ref_fasta,
               ref_fai  = ref_fai,
               ref_amb = ref_amb,
               ref_ann = ref_ann,
               ref_bwt = ref_bwt,
               ref_pac = ref_pac,
               ref_sa = ref_sa,
           sample_name = sample_name,
           library = flowcell[0],
           flowcell = flowcell[1],
           lane = flowcell[2],
           fastq_1 = flowcell[3],
           fastq_2 = flowcell[4]
        }
    }
            call MergeBams as merge {
              input: 
                 sample = sample_name,
                 inBams = align.sorted_bam
        }
            call Artifacts {
       input:
          in_bam = merge.out_bam,
          in_bam_index = merge.out_bam_index,
          sample_name = sample_name,
          reference = ref_fasta
        }

}

task AlignAndSort {
  input {
    File ref_fasta
    File ref_fai
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    
    String sample_name
    String flowcell
    String lane
    String library
    File fastq_1
    File fastq_2

  }
  output {
    File sorted_bam = "~{sample_name}.~{flowcell}.~{lane}.~{library}.sort.bam"
  }
   runtime {
   docker: "gdaly9000/mitochondrial"
   cpu: 32
   memory: "64 GB"
   disks: "local-disk 200 SSD"
   preemptible: 1
 }
  
  command {
    set -e
    samtools --version
    echo "BWA MD5Sum" $(md5sum $(which bwa))

    bwa mem -t 16 -M -R "@RG\\tID:~{flowcell}.~{lane}\\tPU:~{flowcell}.~{lane}.~{sample_name}\\tSM:~{sample_name}\\tPL:ILLUMINA\\tLB:~{library}" \
  ~{ref_fasta} ~{fastq_1} ~{fastq_2}  | \
  samtools sort -n --threads 16 -O BAM -o name.bam -
  samtools fixmate --threads 4 -c -u -O bam name.bam - | \
  samtools sort --threads 28 -O BAM -o ~{sample_name}.~{flowcell}.~{lane}.~{library}.sort.bam - 
  
  }
}

task MergeBams {
   input{
      String sample
      Array[File] inBams
   }
   output{
   File out_bam = "~{sample}.bam"
   File out_bam_index = "~{sample}.bam.bai"
   File bam_stats = "~{sample}.stats.txt"
   Int raw_read_count = read_int("raw_reads.txt")
   }
   runtime {
   docker: "gdaly9000/mitochondrial"
   cpu: 4
   memory: "8 GB"
   disks: "local-disk 200 SSD"
   preemptible: 1
 }
  
  command {
    set -e
  samtools merge --threads 16 -O BAM "~{sample}.bam" ~{sep=" " inBams}
  samtools index "~{sample}.bam"
  samtools flagstat "~{sample}.bam" > "~{sample}.stats.txt"
  # get file with only raw number of reads
  samtools flagstat "~{sample}.bam" | head -1 | cut -d " " -f 1 > raw_reads.txt
  
  }
}

task Artifacts {
   input{
     String sample_name
     File in_bam
     File in_bam_index
     File reference
   }
   output{
     Array[File] out_artifacts = glob("~{sample_name}.artifacts*")
   }
   runtime {
      docker: "docker.io/broadinstitute/gatk:latest"
      cpu: 4
      memory: "8 GB"
      disks: "local-disk 200 SSD"
      preemptible: 1
   }
   command {
     set -e
     gatk --version

     gatk CollectSequencingArtifactMetrics --java-options "-Xmx8g" --INPUT ~{in_bam} \
        --OUTPUT "~{sample_name}.artifacts" \
        --REFERENCE_SEQUENCE ~{reference} --INCLUDE_DUPLICATES true
        
    gatk ConvertSequencingArtifactToOxoG --INPUT_BASE "~{sample_name}.artifacts"
   }
   
}
