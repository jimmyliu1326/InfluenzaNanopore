# check if file is gz compressed
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

# define input for centrifuge
def centrifuge_input(wildcards):
  sample=wildcards["sample"]
  path=samples_meta.Path[wildcards.sample]
  file=glob.glob(path+"/*.fastq*")[0]
  if config["trim"] == 0:
    if is_gz_file(file):
      return os.path.join(sample, sample+".fastq.gz")
    else:
      return os.path.join(sample, sample+".fastq")
  else:
    return os.path.join(sample, "porechop", sample+"_trimmed.fastq")