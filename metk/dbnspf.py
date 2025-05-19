import os 
import pandas as pd
from tqdm import tqdm 
import logging
import subprocess

logger = logging.getLogger()

def dtype_(x):
    try:
        return float(x)
    except:
        return 0

def get_info(line):
    if line == '.':
        return {}
    else:
        fd = {}
        for i in line.split(';'):
            try:
                value = i.split('=')[1].split(',')[0]
                fd[i.split('=')[0]] =  dtype_(value)
            except:
                pass
        
        return fd
    

class DbnSPF():
    def __init__(self, jar, fields=[]):
        self.jar = jar
        if not fields:
            self.default_fields = [
                'chrom', 
                'start', 
                'ref', 
                'alt', 
                'ID', 
                'dbNSFP_BayesDel_addAF_rankscore', 
                'dbNSFP_BayesDel_noAF_rankscore', 
                'dbNSFP_CADD_raw_rankscore', 
                'dbNSFP_ClinPred_rankscore', 
                'dbNSFP_DANN_rankscore', 
                'dbNSFP_Eigen_PC_raw_coding_rankscore',
                'dbNSFP_Eigen_raw_coding_rankscore',
                'dbNSFP_ExAC_AF',
                'dbNSFP_ExAC_nonTCGA_AF',
                'dbNSFP_FATHMM_converted_rankscore',
                'dbNSFP_GERP___NR',
                'dbNSFP_GERP___RS',
                'dbNSFP_GERP___RS_rankscore',
                'dbNSFP_GM12878_fitCons_score',
                'dbNSFP_GenoCanyon_score',
                'dbNSFP_H1_hESC_fitCons_score',
                'dbNSFP_M_CAP_score',
                'dbNSFP_MVP_rankscore',
                'dbNSFP_MetaLR_rankscore',
                'dbNSFP_MetaSVM_rankscore',
                'dbNSFP_MutationTaster_converted_rankscore',
                'dbNSFP_PROVEAN_converted_rankscore',
                'dbNSFP_REVEL_rankscore',
                'dbNSFP_SIFT4G_converted_rankscore',
                'dbNSFP_SiPhy_29way_logOdds_rankscore',
                'dbNSFP_VEST4_rankscore',
                'dbNSFP_bStatistic_converted_rankscore',
                'dbNSFP_fathmm_MKL_coding_rankscore',
                'dbNSFP_fathmm_XF_coding_rankscore',
                'dbNSFP_fathmm_XF_coding_score',
                'dbNSFP_gnomAD_exomes_AF',
                'dbNSFP_integrated_fitCons_rankscore',
                'dbNSFP_phastCons100way_vertebrate',
                'dbNSFP_phastCons30way_mammalian',
                'dbNSFP_DEOGEN2_rankscore',
                'dbNSFP_LRT_converted_rankscore',
                'dbNSFP_MutPred_score',
                'dbNSFP_MutationAssessor_rankscore',
                'dbNSFP_Polyphen2_HDIV_rankscore',
                'dbNSFP_Polyphen2_HVAR_rankscore',
            ]
        else:
            self.default_fields = fields
        
    def format_call(self, annfile):
        scores = []
        for i in tqdm(open(annfile)):
            if "#" == i[0]:
                continue

            chrom, start, _id, ref, alt, _, _, info = i.strip().split()
            score = get_info(info)
            score.update({"chrom": chrom, "start": start, "ref": ref, "alt": alt, "ID": "VAR_{}".format(_id) })

            scores.append(score)
            
        table = pd.DataFrame(scores).fillna(0)
        
        return table

    def match_ids(self, table, vcf):
        raw_table = pd.read_csv(vcf, sep='\t')
        raw_table["ID"] = ["VAR_{}".format(i) for i in raw_table['ID']]
        
        ids = raw_table.loc[:, ['ID']]
        
        table = pd.merge(ids, table, on='ID', how='left').fillna(0)

        # Check if the default fields exists. If not then add 0's. 
        for field in self.default_fields:
            try:
                table[field]
            except Exception as Inst:
                print('This field {} does not match your input data. Adding Zeros to it {}'.format(field, Inst))
                table[field] = 0

        return table[self.default_fields]
    
    def score(self, vcf='', db='', outpath=''):
        self.db = db
        self.outpath = outpath
        
        self.fields = 'hg18_chr' # Added as a placeholder for not getting this field. It can be any unused field reported by snpSift

        if not self.outpath:
            logger.info("Sorting input file ...")
            cmd = 'sort -k1,1 -k2,2n {vcf} > {vcf}.sorted'.format(vcf=vcf)
            os.system(cmd)
            
            logger.info("Matching input file to DBSNPF {} ...".format(self.db))
            cmd = "java -jar {jar} dbnsfp -v -n -f {fields} -db {db} {vcf}.sorted > {vcf}.dbnspf.ann".format(fields=self.fields, jar=self.jar, db=self.db, vcf=vcf)
            logger.info('Running snpsift command: {}'.format(cmd))

            proc = subprocess.Popen(cmd,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        )
            stdout_value, stderr_value = proc.communicate()
            logger.info("dbnSPF log: {}".format(stderr_value))

            # p = subprocess.Popen([cmd], stdout=subprocess.PIPE)
            # logger.info("dbnSPF log: {}".format(p.stdout.read()))
            
            logger.info("Building feature table ...")
            table = self.format_call( "{vcf}.dbnspf.ann".format(vcf=vcf) )
            
            return self.match_ids(table, vcf)
            
        else:
            logger.info("Sorting input file ...")
            cmd = 'sort -k1,1 -k2,2n {vcf} > {vcf}.sorted'.format(vcf=vcf)
            os.system(cmd)
            
            logger.info("Matching input file to DBSNPF {} ...".format(self.db))
            cmd = "java -jar {jar} dbnsfp -v -n -f {fields} -db {db} {vcf}.sorted > {vcf}.dbnspf.ann".format(fields=self.fields, jar=self.jar, db=self.db, vcf=vcf)
            logger.info('Running snpsift command: {}'.format(cmd))

            proc = subprocess.Popen(cmd,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        )
            stdout_value, stderr_value = proc.communicate()
            logger.info("dbnSPF log: {}".format(stderr_value.decode('utf-8')))

            # p = subprocess.Popen([cmd], stdout=subprocess.PIPE)
            # logger.info("dbnSPF log: {}".format(p.stdout.read()))
            
            logger.info("Building feature table ...")
            table = self.format_call( "{out}/{oname}.dbnspf.ann".format(out=self.outpath, oname=vcf.split('/')[-1]) )

            return self.match_ids(table, vcf)