# METK: Mutational Enrichment Toolkit

## Usage
```python
from metk.dataset import read, Dataset
from metk.extractor import FeatureExtractor
import pickle
```

### Preprocess dataset
We need to make sure to have our dataset in the right format. METK uses the convention from cbiportal as the standard way to store the data. If your input dataset is not in the cbio format, you need to have at least the following minimal data schema (columns). 

>**⚠️ Warning:** Your data has to be 0-based following cbioportal convention, if it is not, the results will not be correct. 

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
    mutation_model = 'mutation_model.bin'
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