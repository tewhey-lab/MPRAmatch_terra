# Pipeline for matching barcodes and oligos for MPRA data
# output of the Parse task should be used as the input for the MPRAcount pipeline

workflow MPRAmatch {
  File read_a #R1 fastq
  File read_b #R2 fastq
  File reference_fasta #Oligo sequences with names (can be the oligo order sheet)
  File? attributes #Optional attributes file used for saturation mutagenesis libraries
  Int? barcode_orientation = 2 #2 if you followed the method above, otherwise 1. Default to 2
  Int? thread = 30 #Number of threads to be passed to FLASH2 and MiniMap2. Default to 30
  Int? mem = 30 #Memory to be passed to the sort function. Default to 30G
  Int? read_len = 250 #Length of reads that are being flashed. If mixed lengths use max.
  Int? frag_len = 274 #Fragment length passed to FLASH2
  Int? seq_min = 100 #Minimum acceptable sequence length when separating the barcodes and oligos
  Int? enh_min = 50 #Minimum acceptable length for an oligo
  Int? enh_max = 210 #Maximum acceptable length for an oligo
  Int? bc_len = 20 #Length of barcodes used for project
  String? docker_tag = "latest" #String of the directory relative to the WDL where the other required scripts live
  #String out_directory #String of the directory that all files will be copied to
  String id_out #Project identifier - all files will have this as the prefix for their name
  String? barcode_link = "TCTAGA" #6 base sequence on the barcode end of the link between the barcode and oligo - orientation barcode to oligo
  String? oligo_link = "AGTG" #4 base sequence on the oligo end of the link between the barcode and oligo - orientation barcode to oligo
  String? end_oligo_link = "CGTC" #4 base sequence at the very end of the oligo

  Int disk_pad = 7

  Int flash_disks = ceil(size(read_a, "GB")*1.8 + size(read_b, "GB")*1.8) + disk_pad
  Int pull_disks = ceil(4.25*size(Flash.out, "GB")) + disk_pad
  Int rearr_disks = ceil(1.5*size(Pull_Barcodes.out1, "GB")) + disk_pad
  Int map_disks = ceil(13*size(Rearrange.out, "GB")) + disk_pad
  Int sam2_disks = ceil(1.9*size(MiniMap.out1, "GB")) + disk_pad
  Int sort_disks = ceil(2*size(SAM2MPRA.out, "GB")) + disk_pad
  Int ct_disks = ceil(size(Sort.out, "GB")) + disk_pad
  Int parse_disks = ceil(2*size(Ct_Seq.out, "GB")) + ceil(0.5*disk_pad)
  Int qc_disks = disk_pad
  Int pre_disks = disk_pad

  Int sort_mem = ceil(2*mem)

  call Flash { input:
                  read_a=read_a,
                  read_b=read_b,
                  flash_thread=thread,
                  read_len=read_len,
                  frag_len=frag_len,
                  id_out=id_out,
                  docker_tag=docker_tag,
                  flash_disks=flash_disks
                }
  call Pull_Barcodes { input:
                          #pull = pull,
                          docker_tag=docker_tag,
                          merged_fastq=Flash.out,
                          read_number=barcode_orientation,
                          id_out=id_out,
                          barcode_link=barcode_link,
                          seq_min=seq_min,
                          enh_min=enh_min,
                          enh_max=enh_max,
                          bc_len=bc_len,
                          oligo_link=oligo_link,
                          end_oligo_link=end_oligo_link,
                          pull_disks=pull_disks
                        }
  call Rearrange { input:
                      matched_barcodes=Pull_Barcodes.out1,
                      id_out=id_out,
                      docker_tag=docker_tag,
                      rearr_disks=rearr_disks
                    }
  call MiniMap { input:
                    reference_fasta=reference_fasta,
                    organized_fasta=Rearrange.out,
                    map_thread=thread,
                    id_out=id_out,
                    docker_tag=docker_tag,
                    map_disks=map_disks
                  }
  call SAM2MPRA { input:
                      #sam=sam,
                      docker_tag=docker_tag,
                      sam_file=MiniMap.out1,
                      id_out=id_out,
                      sam2_disks=sam2_disks
                    }
  call Sort { input:
                  MPRA_out=SAM2MPRA.out,
                  mem=mem,
                  id_out=id_out,
                  docker_tag=docker_tag,
                  sort_disks=sort_disks,
                  sort_mem=sort_mem
                }
  call Ct_Seq { input:
                    #count=count,
                    docker_tag=docker_tag,
                    sorted=Sort.out,
                    id_out=id_out,
                    ct_disks=ct_disks
                  }
  if (defined(attributes)) {
    call Parse_sat_mut { input:
                  #  parse=parse,
                    attributes=attributes,
                    docker_tag=docker_tag,
                    counted=Ct_Seq.out,
                    id_out=id_out
                  }
  }
  if (!defined(attributes)) {
    call Parse { input:
                  #  parse=parse,
                    docker_tag=docker_tag,
                    counted=Ct_Seq.out,
                    id_out=id_out,
                    parse_disks=parse_disks
                  }
  }

  call preseq { input:
                 counted=Ct_Seq.out,
                 id_out=id_out,
                 docker_tag=docker_tag,
                 pre_disks=pre_disks
              }
  call qc_plot_t { input:
                    parsed=Parse.out_parsed,
                    hist=Parse.out_hist,
                    preseq_out=preseq.res,
                    preseq_in=preseq.hist,
                    reference_fasta=reference_fasta,
                    docker_tag=docker_tag,
                    id_out=id_out,
                    qc_disks=qc_disks
              }
  #call relocate { input:
  #                  flashed=Flash.out,
  #                  matched=Pull_Barcodes.out1,
  #                  rejected=Pull_Barcodes.out2,
  #                  organized_fasta=Rearrange.out,
  #                  sam_file=MiniMap.out3,
  #                  map_log=MiniMap.out2,
  #                  MPRA_out=SAM2MPRA.out,
  #                  sorted=Sort.out,
  #                  counted=Ct_Seq.out,
  #                  parsed=Parse.out_parsed,
  #                  plothist=Parse.out_hist,
  #                  preseq_hist=preseq.hist,
  #                  preseq_res=preseq.res,
  #                  qc_plot=qc_plot_t.plots,
  #                  out_directory=out_directory
  #                }
  }

task Flash {
  # Flashing raw fastq files together
  File read_a
  File read_b
  Int read_len
  Int frag_len
  Int flash_thread
  Int flash_disks
  String id_out
  String docker_tag

  command {
    flash2 -z -r ${read_len} -f ${frag_len} -s 25 -o ${id_out}.merged -t ${flash_thread} ${read_a} ${read_b}
    }
  output {
    File out="${id_out}.merged.extendedFrags.fastq.gz"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "15G"
    cpu: 32
    disks: "local-disk ${flash_disks} SSD"
    }
  }
task Pull_Barcodes {
  # Pull barcodes from the barcode oligo sequences
  File merged_fastq
  Int read_number
  String docker_tag
  String id_out
  String barcode_link
  String oligo_link
  String end_oligo_link
  Int seq_min
  Int enh_min
  Int enh_max
  Int bc_len
  Int pull_disks
  command <<<
    gzip -dc ${merged_fastq} | perl /scripts/pull_barcodes.pl - ${read_number} ${id_out}.merged ${barcode_link} ${oligo_link} ${end_oligo_link} ${seq_min} ${enh_min} ${enh_max} ${bc_len}
    >>>
  output {
    File out1="${id_out}.merged.match"
    File out2="${id_out}.merged.reject"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${pull_disks} SSD"
    }
  }
task Rearrange {
  # Rearrange the ouptut of the pull_barcodes task to be in a fasta format
  File matched_barcodes
  String id_out
  String docker_tag
  Int rearr_disks
  command <<<
    awk '{print ">"$1"#"$3"\n"$4}' ${matched_barcodes} | gzip - > ${id_out}.merged.match.enh.fa.gz
    >>>
  output {
    File out="${id_out}.merged.match.enh.fa.gz"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${rearr_disks} SSD"
    }
  }
task MiniMap {
  # Map the oligos to the reference to get the oligo names
  File reference_fasta
  File organized_fasta
  Int map_thread
  Int map_disks
  String id_out
  String docker_tag
  command {
    minimap2 --for-only -Y --secondary=no -m 10 -n 1 -t ${map_thread} --end-bonus 12 -O 5 -E 1 -k 10 -2K25m --eqx --cs=short -c -a ${reference_fasta} ${organized_fasta} > ${id_out}.merged.match.enh.sam
    samtools view -S -b ${id_out}.merged.match.enh.sam > ${id_out}.merged.match.enh.bam
    }
  output {
    File out1="${id_out}.merged.match.enh.sam"
    File out3="${id_out}.merged.match.enh.bam"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "200G"
    cpu: 32
    disks: "local-disk ${map_disks} SSD"
    }
  }
task SAM2MPRA {
  # Convert the output of the SAM file into a format closer to the dictionary. Pulls the CIGAR and error information
  File sam_file
  String docker_tag
  String id_out
  Int sam2_disks
  command {
    perl /scripts/SAM2MPRA_cs.pl -C ${sam_file} ${id_out}.merged.match.enh.mapped
    }
  output {
    File out="${id_out}.merged.match.enh.mapped"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "6G"
    cpu: 8
    disks: "local-disk ${sam2_disks} SSD"
    }
  }
task Sort {
  File MPRA_out
  Int mem
  Int sort_disks
  Int sort_mem
  String id_out
  String docker_tag
  command {
    sort -S${mem}G -k2 ${MPRA_out} > ${id_out}.merged.match.enh.mapped.barcode.sort
    }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.sort"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "${sort_mem}G"
    disks: "local-disk ${sort_disks} SSD"
    }
  }
task Ct_Seq {
  # Counts the number of times a barcode-oligo pair occurs
  File sorted
  String docker_tag
  String id_out
  Int ct_disks
  command {
    perl /scripts/Ct_seq.pl ${sorted} 2 4 > ${id_out}.merged.match.enh.mapped.barcode.ct
    }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.ct"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${ct_disks} SSD"
    }
  }
task Parse {
  # Parses the barcode oligo pairs to resolve multimapping barcodes
  File counted
  String docker_tag
  String id_out
  Int parse_disks
  command <<<
    perl /scripts/parse_map.pl ${counted} > ${id_out}.merged.match.enh.mapped.barcode.ct.parsed

    awk '($5 == 0)' ${counted} | awk '{ct[$2]++;cov[$2]+=$4;} END {for(i in ct) print i"\t"ct[i]"\t"cov[i]}' > ${id_out}.merged.match.enh.mapped.barcode.ct.plothist
    >>>
  output {
    File out_parsed="${id_out}.merged.match.enh.mapped.barcode.ct.parsed"
    File out_hist="${id_out}.merged.match.enh.mapped.barcode.ct.plothist"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${parse_disks} SSD"
    }
  }

task Parse_sat_mut {
  # Parses the barcode oligo pairs to resolve multimapping barcodes
  File? attributes
  File counted
  String docker_tag
  String id_out
  Boolean sat_mut=defined(attributes)
  command <<<
    perl /scripts/parse_map.pl -S -A ${attributes} ${counted} > ${id_out}.merged.match.enh.mapped.barcode.ct.parsed
    awk '($5 == 0)' ${counted} | awk '{ct[$2]++;cov[$2]+=$4;} END {for(i in ct) print i"\t"ct[i]"\t"cov[i]}' > ${id_out}.merged.match.enh.mapped.barcode.ct.plothist
    >>>
  output {
    File out_parsed="${id_out}.merged.match.enh.mapped.barcode.ct.parsed"
    File out_hist="${id_out}.merged.match.enh.mapped.barcode.ct.plothist"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk 50 SSD"
    }
  }
task qc_plot_t {
  File parsed
  File hist
  File preseq_out
  File preseq_in
  File reference_fasta
  String docker_tag
  String id_out
  Int qc_disks
  command {
    Rscript /scripts/mapping_QC_plots.R ${parsed} ${hist} ${preseq_out} ${preseq_in} ${reference_fasta} ${id_out}
    }
  output {
    File plots="${id_out}_barcode_qc.pdf"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${qc_disks} SSD"
    }
  }
task preseq {
  # Determine sequencing depth
  File counted
  String id_out
  String docker_tag
  Int pre_disks
  command <<<
    awk '{ct[$4]++}END{for (i in ct)print i "\t" ct[i]}' ${counted} | sort -k1n > ${id_out}.merged.match.enh.mapped.barcode.ct.hist
    preseq lc_extrap -H ${id_out}.merged.match.enh.mapped.barcode.ct.hist -o ${id_out}.merged.match.enh.mapped.barcode.ct.hist.preseq -s 25000000 -n 1000 -e 1000000000
    >>>
  output {
    File hist="${id_out}.merged.match.enh.mapped.barcode.ct.hist"
    File res="${id_out}.merged.match.enh.mapped.barcode.ct.hist.preseq"
   }
  runtime {
    docker: "quay.io/tewhey-lab/mpramatch:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${pre_disks} SSD"
    }
 }
#task relocate{
# File flashed
# File matched
# File rejected
# File organized_fasta
# File sam_file
# File map_log
# File MPRA_out
# File sorted
# File counted
# File parsed
# File plothist
# File preseq_hist
# File preseq_res
# File qc_plot
# String out_directory
# command {
#     mv ${flashed} ${matched} ${rejected} ${organized_fasta} ${sam_file} ${map_log} ${MPRA_out} ${sorted} ${counted} ${parsed} ${plothist} ${preseq_hist} ${preseq_res} ${qc_plot} ${out_directory}
#   }
# }
