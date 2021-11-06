version 1.0

workflow coverageAndInsert {

   input {
    File sample_name
    File bam_file
    File bam_index
    
    Array[File] bed_names
    Array[File] bed_files
    
    File fasta
    File fasta_index
    }
    output {
    Array[File] coverage_stats = process_beds.coverage_stats
    Array[File] coverages = process_beds.coverages
    Array[File] insert_histograms = process_beds.insert_histogram
    Array[File] insert_stats = process_beds.insert_stats

    }

    Array[Int] bed_indices = range(length(bed_names))    
    scatter (bed_index in bed_indices) {

        call Postprocess as process_beds{
           input:
               sample_name = sample_name,
               sample_bam = bam_file,
               sample_index  = bam_index,
               groupName = bed_names[bed_index],
               bed = bed_files[bed_index]
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
    #Array[File] output_tsvs = glob("*/*.tsv")
    File coverage_stats = "coverage-outputs/*.coverage.stats.tsv"
    File coverages = "coverage-outputs/*.coverage.tsv"
    File insert_stats = "insert-outputs/*.insert.stats.tsv"
    File insert_histogram = "insert-outputs/*.insert.hist.tsv"
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
