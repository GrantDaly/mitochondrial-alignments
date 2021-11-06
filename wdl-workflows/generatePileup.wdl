version 1.0

workflow GeneratePileup {

   input {
    File sample_name
    File bam_file
    File bam_index
    
    File fasta
    File fasta_index
    }
    output {
    File raw_vcf = pileup.raw_vcf
    File raw_vcf_index = pileup.raw_vcf_index
    
    }

        call GenerateRawVCF as pileup {
           input:
               sample_name = sample_name,
               sample_bam = bam_file,
               sample_index  = bam_index,
               fasta = fasta,
               fasta_index = fasta_index
        }
    }

task GenerateRawVCF {
   input {
      String sample_name
      File sample_bam
      File sample_index
      File fasta
      File fasta_index
   }
   output {
      File raw_vcf = "~{sample_name}.raw.vcf.gz"
      File raw_vcf_index = "~{sample_name}.raw.vcf.gz.csi"

   }
   runtime {
   docker: "gdaly9000/mitochondrial"
   cpu: 4
   memory: "32 GB"
   disks: "local-disk 200 SSD"
 }
 command {
   bcftools --version
   bcftools mpileup -f ~{fasta} -r chrM -d 10000 -q 20 -Q 20 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR ~{sample_bam}| bcftools call -m -A | bcftools filter -e 'ALT=="."' -O z -o "~{sample_name}.raw.vcf.gz"
   bcftools index "~{sample_name}.raw.vcf.gz"
   
 }
}

 
