import re
import subprocess 
import os 

from rich import print
from collections import Counter
import numpy as np 
import pandas as pd
import logging

logger = logging.getLogger()

class SnpEff():
    
    ''' Extract features from VCF file using SNPEFF '''
    
    def __init__(self, jar='', database=''):
        
        self.ann_fields = [i.strip() for i in 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'.split('|')]
        self.lof_fields = [i.strip() for i in 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'.split('|')]
        self.nmd_fields = [i.strip() for i in 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'.split('|')]

        self.annotation_fields = dict(ANN=self.ann_fields, LOF=self.lof_fields, NMD=self.nmd_fields)
        
        self.bin = jar
        self.database = database
    
    def info2json(self, info):
        fields = {'ANN': [], 'LOF': [], 'NMD': []}
        for i in info.split(';'):
            try:
                field_type, items = i.split('=')
                _items = []
                for item in items.split(','):
                    _items.append({i:j for i,j in zip(self.annotation_fields[field_type], item.split('|')) })

                fields[field_type] = _items
            except:
                pass

        return fields

    def norm_feature(self, feat):
        total = np.sum(list(feat.values()))
        return {i:j/total for i,j in feat.items()}

    def get_conservation_scores(self, mutations):
        econs = []
        for mutation in mutations:
            gene, ref, position, alt = mutation
            prot_name = gene2prot[gene]

            try:
                econs.append(LIST_DC['{}_{}'.format(prot_name, position)])
            except:
                pass
        return econs

    def extract_features(self, mut):
        ann = []
        ann_i = []
        prot = []
        plist = []
        genelist = []

        unk = {}
        for i in mut['ANN']:
            ann.append(i['Annotation'])
            ann_i.append(i['Annotation_Impact'])    

            try:
                protein_info = re.split(r"([0-9]+)", "{}".format(i['HGVS.p']).replace('p.',''))
                prot.append("{}>{}".format(protein_info[0], protein_info[2]))

            except:
                pass

        try:
            lof = float(mut['LOF'][0]['Percent_of_transcripts_affected'].replace(')', ''))
        except:
            lof = 0

        try:
            nmd = float(mut['NMD'][0]['Percent_of_transcripts_affected'].replace(')', ''))
        except:
            nmd = 0

        ann=self.norm_feature(Counter(ann))
        prot=self.norm_feature(Counter(prot))
        ann_i=self.norm_feature(Counter(ann_i))

        out = {"LOF": lof}
        out.update({'NMD':nmd})
        out.update(ann)
        out.update(prot)
        out.update(ann_i)

        return out
    
    def call_snpeff(self, infile='', database=''):
        
        cmd = ['java', '-jar', self.bin, database, infile, '>', '{}.ann'.format(infile)]
        logger.info('Running command: {}'.format(cmd))

        proc = subprocess.Popen(" ".join(cmd),
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        )
        stdout_value, stderr_value = proc.communicate()

        logger.info("snpEff log: {}".format(stderr_value.decode('utf-8')))
        
    def call(self, infile='', database=''):
        self.call_snpeff(infile=infile, database=database)

        return self.load_results(infile=infile, database=database)

    def load_results(self, infile='', database=''):
        features = []
        for ix,i in  enumerate(open('{}.ann'.format(infile))):
            if i[0] == '#':
                continue

            chrom, position, var_id, ref, alt, _, _, info  = i.strip().split('\t')
            infos = self.info2json(info)
            feat = self.extract_features(infos)
            feat.update({'var_id': var_id})

            features.append(feat)

        return pd.DataFrame(features).fillna(0)      