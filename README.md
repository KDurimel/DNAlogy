<h1 align="center">
  <img height="150" src="https://github.com/KDurimel/GSEAtools/blob/master/logo.png" alt="DNAlogy logo"/>
  <br/> 
  DNAlogy
</h1>

# Light bioinformatics tools for gene set statistical analysis and management based on DNA-seq (meta) data

<img src="https://github.com/KDurimel/DNAlogy/blob/master/sep.png" alt="sep"/>

  * [Tool 1 - fastGSEA](#tool-1---fastgsea)
    + [Installation](#installation)
    + [Usage and examples](#usage-and-examples)
    + [Workflow](#workflow)
    + [Methods](#methods)
    + [Contributions](#contributions)
  * [Tool 2 - Javascript API](#tool-2---javascript-api)


## Tool 1 - fastGSEA


[FastGSEA](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA) (fast Gene Set Enrichment Analysis)  performs GO-terms enrichment analysis bewteen two gene sets, based on hypergeometric tests. These gene sets must be provided as text files containing one international databank (ncbi, refseq, etc) gene or proteins identifier per line. **FastGSEA can also be used only as a standalone mapping tool**, using the `--mapOnly` option.

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

### Installation

Only on GNU/Linux based systems. You can install it using `Conda`. If you already have Conda on your system you can apply the following instructions:

**Automatically**

Just go to the FAST_GSEA directory, then run `install.sh`

```bash
# Go to fastGSEA directory
cd DNAlogy/FAST_GSEA

# Run install script
bash install.sh 

# Respond "yes" to the "Do you wish the installer to prepend the Miniconda2 install location to PATH in your /home/username/.bashrc ?" answer.

# Activate the enrivonment :
source activate gsea_env
```
Next, you'll just have to activate the `gsea_env` environment when you want to use fastGSEA.

```bash
source activate gsea_env
```

**Manually**

Use the [packages.yml](https://github.com/KDurimel/DNAlogy/tree/master/FAST_GSEA/packages.yml) to retrieve all Python and R dependencies then install fastGSEA :

```bash
# Go to fastGSEA directory
cd DNAlogy/FAST_GSEA

# Create the environment from the yaml file
conda env create -f packages.yml 

# Activate the enrivonment :
source activate gsea_env

# Add fastGSEA to your environment variables
PATH="src/:${PATH}"
export PATH

# or add an alias for fastGSEA in your bashrc
echo 'alias fastGSEA="python $(pwd)'/src/fastGSEA.py'"' >> ~/.bashrc

# Reload your shell (here, bash) settings 
source ~/.bashrc
```


------

### Usage and examples


**Command line options**

| Option        | Description                                                             | Required |
|---------------|-------------------------------------------------------------------------|----------|
| -ech          | Sample ids file (one id per line)                                       | Yes      |
| -univ         | Universe ids file (one id per line)                                     | Yes      |
| -mappingFile  | idmapping_selected.tab.gz file                                          | Yes      |
| -output       | Output results prefix                                                   | Yes      |
| -obo          | Gene ontology .obo graph file used when "--trim" option activated       | No       |
| -toDB         | databank identifier wanted as output when "--mapOnly"  option activated | No       |
| --fromOtherDB | Activate all ids support (slower)                                       | No       |
| --mapOffline  | Perform "MAP" step offline                                              | No       |
| --trim        | Trim prokaryotic GO-terms                                               | No       |
| --view        | Outputs the enrichment results in a 2D graph                            | No       |
| --mapOnly     | Perform only the "MAP" step and keep its results                        | No       |
| --keepTmp     | Keep temporary files folder                                             | No       |




**Usage - Where can I find the `-mappingFile` and `-obo` files?**

`-mappingFile` is used for offline ids mappping and can be found on the [Uniprot FTP](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz) ( > 6gb file)

`-obo` is the gene ontology graph used for GO-terms checking and trimming. It is availaible on the [Open Biological and Biomedical Ontology (OBO)](purl.obolibrary.org/obo/go/releases/2018-06-01/subsets/gosubset_prok.obo) ( > 60 mo file). There is daily releases, so you can download the latest ones [here](purl.obolibrary.org/obo/go/releases)


> **Important note :** be careful when trimming non prokaryotic GO-terms, 'gosubset_prok' terms are not maitained since 2018/06 because some of them muight be irrelevant. More information [here](https://github.com/geneontology/go-ontology/pull/16255) and [here](https://github.com/geneontology/go-ontology/issues/16077).



**Usage - Id mapping : map any ids to UniRef100 ids**

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


**Example - Id mapping : map any ids to UniRef100 ids**

This command line requires two outputs, if you want to perform id mapping in only one file, just provide it twice as `-ech` and `-univ`.
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
fastGSEA -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz --mapOnly -toDB UniRef100 -output maybe/here
```


**Examples - Gene set enrichment analysis : find which GO-terms from a gene set are overrepresented**

With most steps offline (faster, the better updated your -mappingFile and/or -obo are, the better the results will be) :
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
fastGSEA -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz  --mapOffline -output maybe/here
```

Requesting NCBI and Uniprot APIs (most reliable, but slower):
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
fastGSEA -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz  -output maybe/here
```

...plus trimming obsolete and non Prokaryotic GO-terms (up to date obo file `gosubset_prok.obo` needed):
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
fastGSEA -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz -obo gosubset_prok.obo -output somewhere --trim
```

...plus generating a chart for enriched GO-terms:
```bash
# -ech: gene set sample ; -univ: gene set universe ; other args? please read the docs :)
fastGSEA -ech ech.txt  -univ univ.txt  -mappingFile idmapping_very_light.gz -obo gosubset_prok.obo -output somewhere --trim --view
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

## Tool 2 - Javascript API

~ Q4 2018 (W.I.P)
