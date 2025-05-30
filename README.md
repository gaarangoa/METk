# METK: Mutational Enrichment Toolkit

## Usage
```python
from metk.dataset import read, Dataset
from metk.extractor import FeatureExtractor
import pickle
```

### Preprocess dataset
We need to make sure to have our dataset in the right format. METK uses the convention from cbiportal as the standard way to store the data. If your input dataset is not in the cbio format, you need to have at least the following minimal data schema (columns). 

>**⚠️ IMPORTANT:** Make sure your data loci is correct. To do so see **validate my input data** section.

```python 
# Define table schema
schema = {
    "Chromosome": 'Chromosome',
    "Start_Position": 'Start_Position',
    "End_Position": 'End_Position',
    "Reference_Allele": 'Reference_Allele',
    "Tumor_Allele_1": 'Tumor_Seq_Allele1',
    "Tumor_Allele_2": 'Tumor_Seq_Allele2',
    "Variant_Type": 'Variant_Type',
    "Gene_Name": 'Hugo_Symbol',
    "Sample_ID": 'Tumor_Sample_Barcode'
}

# Load data into a dataframe
table = read(
    '../data/brca/brca_tcga.mutations.txt', 
    sep='\t', 
    schema=schema, 
    reference_genome='GRCh37',
    added_fields = []
)

# Remove non SNV/INDEL mutations
table = table[
    ((table.Variant_Type == 'SNV') | (table.Variant_Type == 'SNP') | (table.Variant_Type == 'DEL') | (table.Variant_Type == 'INS'))
].reset_index(drop=True)

table = Dataset(table)
```

### Extract features with METK
Now that  we have cleaned and formatted the data we are ready to extract all the features available from METK.

```python
mutation_enrichment = FeatureExtractor(
    reference_genome = 'GRCh37',
    db = '/mnt/METk/',
    identifier = 'table_unique_id_',
    variant_types = 'INS,DEL,SNP',
    output_path = '../data/brca/',
    mutation_model = 'dgv2.cbioportal.128.e500.bin',
    run_snpeff=False, # set to True if required
    run_dbnsfp=False, # set to True if required
)
mutation_enrichment.extract_features(table)
```

### Post Process 
Features are calculated by each entry in your data (each mutation) and are stored as a pickle file. You can access the 
data by loading this file. There are three objects: 

* vectors: it contains the extracted features for each mutation
* feature_names: The name of the extracted features: 
    * dg_*: deepgesture embeddings
    * ge_*: Embeddings from the gene name
    * dbNSFP: Scores from the dbNSFP database
    * other: SnpEff extracted scores
* table: This is your input table 

**Recommended** you can create a master table with all the vectors and input information as follows: 
```python
import pickle
[vectors, feature_names, metadata] = pickle.load(open('../data/brca/mutation_features.pk', 'rb'))

vectors = pd.DataFrame(vectors, columns = feature_names)
data = pd.concat([metadata, vectors], axis=1)
```

### Patient level embeddings
Patient level embeddings refer to aggregating the extracted features to the patients instead of individual mutations. This is practical if you are 
looking for patterns observed in a patient. 

This piece of code will compute the average over of the features for each patient in your data. For the patient embeddings we will ad a p_* to distinguish from the mutation embeddings. 

```python
patient_embeddings = data.groupby(['Sample_ID']).mean().reset_index()[['Sample_ID'] + feature_names]
patient_embeddings.columns = ['Sample_ID'] + ['p_dg_{}'.format(i) for i in feature_names]
```

You can add the patient level embeddings to the master table (if needed) by merging the two tables.

```python
data = pd.merge(data, patient_embeddings, on='Sample_ID', how='left')
```

Finally, once you have this master table, take a look at a couple of patients and verify that you are getting the expected values. Be careful with merging tables, make sure your sample_ids are unique and represent a patient. For instance, if your dataset contains patients with different visits (DNA extractions) you may have several samples, therefore using the above code would mix all the samples in a patient. You need to consider these scenarios when doing the analysis at the patient level. 

## How to validate my input data
METk generates an auxiliary file that has the extracted context from the human genome, the reference allele should 
match the centered allele on the extracted sequence with left and right context. If there is a mismatch, it means that your chromosome positions are shifted and need to be corrected. Here is a code snippet for doing that:

```python
context_sequences = pd.read_csv('../data/brca/metk/metadata.txt', sep='\t')
context = context_sequences[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele_1', 'Tumor_Allele_2', 'context']]
context['context'] = [f"{i[:10]} {i[10:-10]} {i[-10:]}" for i in context.context]
context.head()
```

| Chromosome | Start_Position | Reference_Allele | Tumor_Allele_1 | Tumor_Allele_2 |     left      | center |     right     |
|------------|----------------|------------------|----------------|----------------|---------------|--------|---------------|
| 3          | 132047117      | C                | C              | T              | TCAGGTGTTT    | C      | GGCATGGAGA    |
| X          | 110951473      | T                | T              | C              | ACCCCCACCC    | T      | GTACAAAATG    |
| 1          | 104166496      | T                | T              | C              | AGGATGTTAA    | T      | GATTGGGTTG    |
| 9          | 33135286       | T                | T              | C              | GGTTGCGGAA    | T      | GGAATGATGA    |
| 5          | 40841761       | A                | A              | C              | TTGCGGCTTA    | A      | GGCATGGTAA    |
| 19         | 48821779       | G                | G              | A              | TCTGGGCTGC    | G      | CTGATCTGCA    |
| 16         | 68862076       | G                | G              | T              | CACCATCCCA    | G      | TTCTGATTCT    |

## Setup
### Environment 
METK uses a docker image with the relevant libraries installed. 
>**⚠️ Info:** We strongly recommend using docker (or singulairty) for extracting embeddings from METk.

### Start jupyter notebook 
```bash 
docker run -it --rm \
    -p 8888:8888 \
    -v /path/to/metk/:/METk/ \
    -v /path/to/data/:/data/ \
    gaarangoa/chip_classifier:version-4.0.0 \
    jupyter notebook --NotebookApp.default_url=/lab/ --ip=0.0.0.0 --port=8888 --allow-root
```

#### inside the container install METk package
```bash
pip install git+https://github.com/gaarangoa/METk.git
```

### Download METk embeddings
Download the METk models and add them to the metk data path:
* <a href="https://zenodo.org/records/15484265/files/dgv2.cbioportal.128.e500.bin?download=1&preview=1"> METK mutational embeddings model </a>
* <a href="https://zenodo.org/records/15484265/files/genes.tcga.model.d8.tsv?download=1&preview=1"> METK gene embeddings model </a>

### Human reference genome
Download the human reference genomes version and add them to the metk data path: 
* <a href="https://zenodo.org/records/15484265/files/reference_genome.tar.gz?download=1&preview=1"> reference genomes </a>

### SnpEff & SnpSift (optional)
You can skip this step if you are interested on getting only embeddings.

METk also extracts scores functional scores from snpEff (v4.3) and snpSift (v5.0).
You need to download both versions from the <a href="https://pcingola.github.io/SnpEff/snpsift/introduction/#download-and-install" open=blank_> official website </a> with the following databases and add them to the metk data path: 
* <a href="https://zenodo.org/records/4323592#.X9zPhNhKiHs">dbNSFP4.1a or dbNSFP4.1c</a> (for snpEff v5.0)
* snpEff_v4_3_GRCh37.75 (for snpEff v4.3)
* snpEff_v4_3_GRCh38.86 (for snpEff v4.3)

To ignore SnpEff/SnpSift set ```run_snpeff=False``` and ```run_dbnsfp=False```.

### METK database structure
Finally, the metk data directory should have the following structure: 
```bash
📁 METk/
|-- 📁 starspace
|   -- genes.tcga.model.d8.tsv
|   -- dgv2.cbioportal.128.e500.bin
|-- 📁 reference_genome
|   -- GRCh37.fa  
|   -- GRCh37.fa.fai  
|   -- GRCh38.d1.vd1.fa  
|   -- GRCh38.d1.vd1.fa.fai
|-- 📁 snpEff (optional)
|   -- 📁 snpEff_5_0
|      -- 📁 db
|         -- 📁 GRCh37
|            -- 📁 dbNSFP
|         -- 📁 GRCh38
|            -- 📁 dbNSFP
|   -- 📁 snpEff_4_3
|      -- 📁 snpEff
|         -- 📁 data
|            -- 📁 GRCh37.75
|            -- 📁 GRCh38.86
```

## Licensing
This software includes third-party components that are governed by separate license terms:

>**⚠️ IMPORTANT:** Users are responsible for reviewing and complying with the licenses of all third-party components used by METk (<a href="https://pcingola.github.io/SnpEff/snpsift/introduction/#download-and-install" open=blank_>snpEff</a>, <a href="https://sites.google.com/site/jpopgen/dbNSFP">Legacy dbNSFP v4.1 (see 1)</a> <a href='https://www.dbnsfp.org/license'>(and 2)</a>). These tools are not distributed with this software and must be installed separately. Their licenses apply independently.
