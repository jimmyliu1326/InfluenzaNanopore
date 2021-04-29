rule read_count:
  input:
    combined_fastq=centrifuge_input,
    segment_fastq=lambda wildcards: expand(wildcards.sample+"/binned_fastq/{segment}.fastq", segment=segment_l)
  output:
    count_summary=temp("{sample}/read_count_summary.csv")
  threads: 4
  shell:
    """
    expr $(wc -l {input.combined_fastq} | cut -f1 -d' ') / 4 > {output.count_summary}
    for fastq in {input.segment_fastq}; do
      expr $(wc -l $fastq | cut -f1 -d' ') / 4  >> {output.count_summary}
    done
    """

rule read_summary:
  input: expand("{sample}/read_count_summary.csv", sample = samples_meta.Sample)
  output: "summary_statistics.csv"
  threads: 1
  params:
    segment=segment_l
  shell:
    """
    echo "sample,total_reads,classified_reads,percent_classified,$(echo {params.segment} | sed 's/ /,/g')" > {output}
    for summary in {input}; do
      sample=$(dirname $summary)
      segment_reads=$(tail -n +2 $summary | paste -s -d',')
      total=$(head -n1 $summary)
      classified=$(tail -n +2 $summary | paste -sd+ | bc)
      percent=$(echo "scale=2;$classified/$total*100" | bc)
      echo "${{sample}},${{total}},${{classified}},${{percent}},${{segment_reads}}" >> {output}
    done
    """

rule read_length_viz:
  input: 
    readlength_res=expand("{sample}/binned_fastq/{segment}.read.lengths", sample=samples_meta.Sample, segment=segment_l),
    summary="summary_statistics.csv"
  output: "InfA_analysis_viz.html"
  params:
    outdir=config["outdir"]
  threads: 4
  script: config["pipeline_dir"]+"/src/InfA_analysis_viz.Rmd"

