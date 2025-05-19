from judge.utils import prog
from rich import print
import os
import pandas as pd
import numpy as np 
import pickle
from collections import Counter

class Dataset(pd.DataFrame):
    def extract(self, ):
        self['Reference_Genome'] = self._reference_genome
        self['table_unique_id_'] = ["uid_{}".format(i) for i in self.index]

        table = self[self.fields.split('\t')]
        table.columns = self.column_names.split('\t')

        duplicated_names = Counter(list(self) + self.column_names.split('\t'))
        duplicated_names = [k for k,v in duplicated_names.items() if v > 1]

        self.drop([ i for i in duplicated_names if i != 'table_unique_id_' ], axis=1, inplace=True)
        
        full_table = pd.merge(table, self, on='table_unique_id_')
        print('Total entries in input table: {}'.format(full_table.shape))

        return full_table
    
    def set_fields(self, ):
        self.fields = "\t".join([self._chromosome, self._start, self._end, self._reference_allele, self._tumor_allele_1, self._tumor_allele_2, 'Reference_Genome', self._variant_type, self._gene_name, self._sample_id, 'table_unique_id_'])
        self.column_names = "\t".join(['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Allele_1', 'Tumor_Allele_2', 'Reference_Genome', 'Variant_Type', 'Hugo_Symbol', 'Sample_ID', 'table_unique_id_'])
    
    def add_column(self, key, new_name=''):
        self.fields += "\t{}".format(key)
        if new_name:
            self.column_names += "\t{}".format(new_name)
        else:
            self.column_names += "\t{}".format(key)
    
    def to_vcf(self, outfile, var_id=''):
        fname = "{}".format(outfile)
        fo = open(fname, 'w')
        fo.write('#CHROM\tPOS\tID\tREF\tALT\n')
        
        for ix, i in self.iterrows():
            _id = ix
            if var_id:
                _id = i[var_id]
            tumor_allele = i.Tumor_Allele_1
            if i.Reference_Allele == i.Tumor_Allele_1:
                tumor_allele = i.Tumor_Allele_2
            
            fo.write("{}\t{}\t{}\t{}\t{}\n".format(i.Chromosome, i.Start_Position, _id, i.Reference_Allele, tumor_allele).replace('-', '.') )


class CustomDataset(Dataset):
    def __init__(self, schema, reference_genome, *args, **kwargs):
        super(Dataset, self).__init__(*args, **kwargs)
        
        self._reference_genome = reference_genome
        self._chromosome = schema["Chromosome"]
        self._start = schema["Start_Position"]
        self._end = schema["End_Position"]
        self._reference_allele = schema["Reference_Allele"]
        self._tumor_allele_1 = schema["Tumor_Allele_1"]
        self._tumor_allele_2 = schema["Tumor_Allele_2"]
        self._variant_type = schema["Variant_Type"]
        self._gene_name = schema["Gene_Name"]
        self._sample_id = schema["Sample_ID"]

def read(infile, schema, reference_genome, added_fields=[], *args, **kwargs):
    
    dataset = CustomDataset(data=pd.read_csv(infile, *args, **kwargs, comment='#', low_memory=False), schema=schema, reference_genome=reference_genome)
    dataset.set_fields()
    
    for i in added_fields:
        dataset.add_column(i, i)
    
    table = dataset.extract()
    
    return Dataset(table)