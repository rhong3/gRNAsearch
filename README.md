# **Search for guide RNA sites in PCR tags**
### Usage
 - Install the package version through pip `pip install gRNAsearch `
 - `from gRNAsearch.search import search`
 - `summary, frames = search(syn="../seq.fa", tag="../tag.csv", save=True, output_dir="~", output_name="gRNA")`
### Parameters 
`search(syn, tag, save=True, output_dir='~', output_name="gRNA"):`
 - syn: full path to sequence file in '.fa' format
 - tag: PCR tag file
 - save: save the outputs in csv (default=True)
 - output_dir: output directory (full path) to save the files (default='~')
 - output_name: name of the output file (default="gRNA")
### Outputs
 - summary: the search results data frame
 - frames: the intermediate data frame with WT sequences
 