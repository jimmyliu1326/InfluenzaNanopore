# define input for centrifuge
def centrifuge_input(wildcards):
  sample=wildcards["sample"]
  if config["trim"] == 0:
    return os.path.join(sample, sample+".fastq")
  else:
    return os.path.join(sample, "porechop", sample+"_trimmed.fastq")