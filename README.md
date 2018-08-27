<h1 align="center">
  <img height="150" src="https://github.com/KDurimel/GSEAtools/blob/master/logo.png" alt="DNAlogy logo"/>
  <br/> 
  DNAlogy
</h1>

# Light bioinformatics tools for gene set statistical analysis and management based on DNA-seq (meta) data

<img src="https://github.com/KDurimel/DNAlogy/blob/master/sep.png" alt="sep"/>

## Tool 1 - fastGSEA (fast Gene Set Enrichment Analysis) 


[FastGSEA](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA) performs GO-terms enrichment analysis bewteen two gene sets, based on hypergeometric tests. These gene sets must be provided as text files containing one international databank (ncbi, refseq, etc) gene or proteins identifier per line. **FastGSEA can also be used only as a standalone mapping tool**, using the `--mapOnly` option.

In every cases of use, supported ids are:

1. UniProtKB-AC 
2. RefSeq 
3. UniProtKB-ID 
4. GeneID (EntrezGene) 
5. GI 
6. GO (as output only!) 
7. UniRef100 
8. UniRef90 
9. UniRef50 
10. UniParc 
11. UniGene 
12. EMBL 
13. EMBL-CDS 
14. Ensembl 
15. Ensembl_TRS 
16. Ensembl_PRO

------

### Installation (GNU/Linux)

You can install it using `Conda`.

For a full installation (all python and R environment), use the [packages_full.yml](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA/packages_full.yml) file:

```bash
# Go to fastGSEA directory
cd FAST_GSEA

# Create the environment from the yaml file
conda env create -f packages_full.yml 

# Activate the enrivonment :
source activate gsea_env # or the first line you just edited

# Add FastGSEA to your environment variables
PATH="src/:${PATH}"
export PATH

# Reload your shell (here, bash) settings 
~/.bashrc
```


For a light installation (python dependencies only) use the [packages.yml](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA/packages.yml) file:

```bash
# Go to fastGSEA directory
cd FAST_GSEA

# Create the environment from the "light" yaml file
conda env create -f packages.yml 

# Activate the enrivonment :
source activate gsea_env # or the first line you just edited

# Add FastGSEA to your environment variables
PATH="src/:${PATH}"
export PATH

# Reload your shell (here, bash) settings 
~/.bashrc
```

------

### Usage and examples

FastGSEA takes two input files (one for sample, second one for universe). They have to be text files containing one international databank (ncbi, refseq, etc) **supported** identifiers (listed above) per line, for example:

```bash
O55719
Q6GZM9
Q6GZM8
NP_302218.1
WP_008262748.1
WP_011437797.1
NP_149806.1
NP_854636.1
WP_003877490.1
Q6GZM7
P0C9F0
P0C9F1
P0C9F2
P0C9E9
Q65209
P0C9F4
P0C9F5
P0C9F6
```
Dummy data and its results can be found in the [example](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA/examples) directory.

**Id mapping : map any ids to UniRef100 ids**

```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
python fastGSEA.py -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz --mapOnly -toDB UniRef100 -output maybe/here
```

**Gene set enrichment analysis : find which GO-terms from a gene set are overrepresented**

With most steps offline (faster, the better updated your -mappingFile and/or -obo are, the better the results will be) :
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
python fastGSEA.py -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz  --mapOffline -output maybe/here
```

Requesting NCBI and Uniprot APIs (most reliable, but slower):
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
python fastGSEA.py -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz  -output maybe/here
```

...plus trimming obsolete and non Prokaryotic GO-terms (up to date obo file `gosubset_prok.obo` needed):
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
python fastGSEA.py -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz -obo gosubset_prok.obo -output somewhere --trim
```

...plus generating a chart for enriched GO-terms:
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
python fastGSEA.py -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz -obo gosubset_prok.obo -output somewhere --trim --view
```
------

All these options can be combined to use FastGSEA as you like. [Dummy data](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA/examples/input_data) can be used to try IT. Examples of results from this data are also provided [here](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA/examples/results)

------

### Workflow
As said previously, workflow can be stopped at each step, the last 3 parts of the workflow are optional and and behave as you set it up for.
<img src="https://github.com/KDurimel/DNAlogy/blob/master/FAST_GSEA/doc/workflow.png" alt="fastGSEA workflow"/>

------

### Methods
FastGSEA comes with several methods that you can manipulate to make it behave as you like. Fore more details, please read the [technical documentation](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA/doc).

------
### Contributions
Opened to all suggestions and pull requests.

<img src="https://github.com/KDurimel/DNAlogy/blob/master/sep.png" alt="sep"/>

## Tool 2 - Javascript API (w.i.p)

~ Q4 2018
