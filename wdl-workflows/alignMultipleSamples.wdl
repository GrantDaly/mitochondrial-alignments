version 1.0

struct References {
    File ref_fasta
    File ref_fai
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
}

struct Flowcell {
  String flowcell_name
  String flowcell_lane
  Array[String] filenames
}

struct Sample {
   String sample_name
   String library_name
   Array[Flowcell] flowcells
}

workflow alignSamples {

   input {
    References references
    Array[Sample] samples
    }
    output {
    Array[File] out_bam = merge.out_bam
    Array[File] out_bam_index = merge.out_bam_index
    Array[Array[File]] out_artifacts = Artifacts.out_artifacts
    }
    
    scatter (sample in samples) {
    scatter (flowcell in sample.flowcells) {
        call AlignAndSort as align{
           input:
               ref_fasta = references.ref_fasta,
               ref_fai  = references.ref_fai,
               ref_amb = references.ref_amb,
               ref_ann = references.ref_ann,
               ref_bwt = references.ref_bwt,
               ref_pac = references.ref_pac,
               ref_sa = references.ref_sa,
           sample_name = sample.sample_name,
           library = sample.library_name,
           flowcell = flowcell.flowcell_name,
           lane = flowcell.flowcell_lane,
           fastq_1 = flowcell.filenames[0],
           fastq_2 = flowcell.filenames[1]
        }
    }
            call MergeBams as merge {
              input: 
                 sample = sample.sample_name,
                 inBams = align.sorted_bam
        }
            call Artifacts {
       input:
          in_bam = merge.out_bam,
          in_bam_index = merge.out_bam_index,
          sample_name = sample.sample_name,
          reference = references.ref_fasta
        }
    }
    #Array[File] merged_bams = merge.out_bam
    #Array[File] merged_bams_index = merge.out_bam_index

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
   }
   runtime {
   docker: "gdaly9000/mitochondrial"
   cpu: 4
   memory: "8 GB"
 }
  
  command {
    set -e
  samtools merge --threads 16 -O BAM "~{sample}.bam" ~{sep=" " inBams}
  samtools index "~{sample}.bam"
  
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
