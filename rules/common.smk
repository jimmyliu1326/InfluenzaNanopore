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

# get segment name
def get_segment_name(segment_str):
  segment_names = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
  segment_idx = int(segment_str.replace("segment_", "")) - 1
  segment_name = segment_names[segment_idx]
  return(segment_name)
  
