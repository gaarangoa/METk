import pickle
import numpy as np
import sys
from tqdm import tqdm
import argparse

from .embeddings import Embedder 
from .snpeff import SnpEff 
from .dbnspf import DbnSPF

import logging

logger = logging.getLogger()

snpeff_items = [
    'LOF',
    'NMD',
    'missense_variant',
    'downstream_gene_variant',
    'MODERATE',
    'MODIFIER',
    'stop_gained',
    'upstream_gene_variant',
    'intron_variant',
    'non_coding_transcript_exon_variant',
    'HIGH',
    '5_prime_UTR_variant',
    'protein_protein_contact',
    'stop_gained&splice_region_variant',
    'structural_interaction_variant',
    'sequence_feature',
    'LOW',
    'missense_variant&splice_region_variant',
    'splice_donor_variant&intron_variant',
    'splice_acceptor_variant&intron_variant',
    '3_prime_UTR_variant',
    'splice_region_variant&non_coding_transcript_exon_variant',
    '5_prime_UTR_premature_start_codon_gain_variant',
    'splice_region_variant&synonymous_variant',
    'splice_region_variant',
    'stop_lost',
    'stop_lost&splice_region_variant',
    'synonymous_variant',
    'splice_region_variant&intron_variant',
    'TF_binding_site_variant',
    'start_lost',
    'stop_retained_variant',
    'start_lost&splice_region_variant',
    'initiator_codon_variant',
    'intragenic_variant',
 ]

snpsift_index = [
    'dbNSFP_BayesDel_addAF_rankscore', 
    'dbNSFP_BayesDel_noAF_rankscore', 
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
    'dbNSFP_H1_hESC_fitCons_score',
    'dbNSFP_M_CAP_score',
    'dbNSFP_MVP_rankscore',
    'dbNSFP_MetaLR_rankscore',
    'dbNSFP_MetaSVM_rankscore',
    'dbNSFP_MutationTaster_converted_rankscore',
    'dbNSFP_PROVEAN_converted_rankscore',
    'dbNSFP_SIFT4G_converted_rankscore',
    'dbNSFP_SiPhy_29way_logOdds_rankscore',
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
    'dbNSFP_MutationAssessor_rankscore',
]



class FeatureExtractor():
    def __init__(self, db, reference_genome, identifier, variant_types, output_path, run_deepgesture=True, run_snpeff=True, run_dbnsfp=True, run_deepgesture_gene=True, **kwargs):
        
        ''' 
        Feature extractor for CHIP classifier:
        Inputs: 
            db: path where the databases are stored ('/data/chip/')
            reference_genome: genome reference to use
            identifier: Unique ID in the table to extract the features. 
            variant_types: names of the variant types to consider. As the names may change depending on the 
                vendor this field is mandatory. For instance, Guardant health reports the variants as 
                Small insertion,Small deletion,Single base substitution. Use the same order if you have different
                names. From TCGA portal it may be: INS,DEL,SNP.
            output_path: Directory where the results will be stored. Use one directory per table. 
        
        Returns: 
            A data frame with the CHIP classifier fields. 
        '''

        self.deepgesture_model = kwargs.get('mutation_model', 'mutation_model.bin')

        self.db = db
        self.reference_genome = reference_genome
        self.variant_types = variant_types
        self.output_path = output_path
        self.identifier = identifier

        self.snpeff_dir = '{}/snpEff/snpEff_4_3/snpEff/snpEff.jar'.format(db)
        self.snpsift_jar = '{}/snpEff/snpEff_5_0/SnpSift.jar'.format(db)
        self.snpsift_db = [
            '{}/snpEff/snpEff_5_0/db/GRCh37/dbNSFP/dbNSFP4.1a.txt.gz'.format(db), 
            '{}/snpEff/snpEff_5_0/db/GRCh38/dbNSFP/dbNSFP4.1a.txt.gz'.format(db)
        ]
        self.gene_embeddings_model = '{}/starspace/genes.tcga.model.d8.tsv'.format(db)

        self.run_deepgesture = run_deepgesture
        self.run_snpeff = run_snpeff
        self.run_dbnsfp = run_dbnsfp
        self.run_deepgesture_gene = run_deepgesture_gene
        
    def compute_deepgesture_embeddings(self, table):
        logger.info('Step 1: Extracting DeepGESTURE embeddings ...')
        
        dg = Embedder(
            path_to_reference_genomes='{}/reference_genome/'.format(self.db),
            path_to_trained_embeddings='{}/starspace/'.format(self.db),
            pretrained_embeddings=self.deepgesture_model
        )
        
        if self.run_deepgesture:
            dg.call(
                table = table, 
                reference_genome = self.reference_genome,
                identifier = self.identifier,
                selected_variant_types=self.variant_types,
                output_path = self.output_path
            )
        
        self.deepgesture_features, self.deepgesture_ids = dg.load_embeddings('{}/embeddings.wv'.format(self.output_path))
        self.deepgesture_names = [ "dg_{}".format(i) for i in range(self.deepgesture_features.shape[1])]

        logger.info('Step 1: Feature dimension: {}'.format(self.deepgesture_features.shape))
        
    def compute_snpeff_features(self, ):
        
        logger.info('Step 2: Extracting SnpEff features ...')
        snpeff = SnpEff(jar=self.snpeff_dir)
        
        database = 'GRCh37.75' if self.reference_genome == 'GRCh37' else 'GRCh38.86'
        
        self.snpeff_features = np.array([])
        if self.run_snpeff:
            snpeff_data = snpeff.call(infile='{}/deepgesture.vcf'.format(self.output_path), database=database) 
            self.snpeff_ids = snpeff_data['var_id']
        
            for key in snpeff_items:
                try:
                    snpeff_data[key]
                except:
                    snpeff_data[key] = 0

            self.snpeff_features = np.array(snpeff_data[snpeff_items], dtype=float)
            self.snpeff_names = snpeff_items

        logger.info('Step 2: Feature dimension: {}'.format(self.snpeff_features.shape))
        
    def compute_dbnsfp_features(self, ):
        logger.info('Step 3: Extracting DbnSPF features ...')
        dbnspf = DbnSPF( jar = self.snpsift_jar, )

        self.snpsift_features = np.array([])
        if self.run_dbnsfp:
            snpsift = dbnspf.score(
                vcf='{}/deepgesture.vcf'.format(self.output_path),
                db=self.snpsift_db[0] if self.reference_genome == 'GRCh37' else self.snpsift_db[1],
                outpath = self.output_path
            )

            self.snpsift_features = np.array(snpsift[snpsift_index], dtype=float)
            self.snpsift_ids = snpsift.iloc[:, :5]
            self.snpsift_names = snpsift_index
            
            logger.info('Step 3: Feature dimension: {}'.format(self.snpsift_features.shape))

    def compute_gene_embeddings(self, table):
        logger.info('Step 4: Extracting Individual Gene features ...')
        stargene = {i.strip().split('\t')[0]:np.array(i.strip().split('\t')[1:], dtype=float) for i in open(self.gene_embeddings_model)}
        
        gene_embeddings = []
        missing_genes = []
        for i in table.Hugo_Symbol:
            try:
                gene_embeddings.append(stargene[i])
            except:
                gene_embeddings.append(np.zeros(8))
                missing_genes.append(i)


        self.gene_features = np.array(gene_embeddings)
        self.gene_feature_names = ["ge_{}".format(i) for i in range(self.gene_features.shape[1])]

        logger.info('Step 4: Feature dimension: {}'.format(self.gene_features.shape))

    def extract_features(self, table): 
        self.compute_deepgesture_embeddings(table)
        self.compute_snpeff_features()
        self.compute_dbnsfp_features()
        self.compute_gene_embeddings(table)

        if self.run_snpeff and self.run_dbnsfp:            
            self.dataset = np.array(np.concatenate([self.deepgesture_features, self.gene_features, self.snpeff_features, self.snpsift_features], axis=1), dtype=float)
        
        elif self.run_snpeff and not self.run_dbnsfp:            
            self.dataset = np.array(np.concatenate([self.deepgesture_features, self.gene_features, self.snpeff_features], axis=1), dtype=float)
        
        elif not self.run_snpeff and self.run_dbnsfp:            
            self.dataset = np.array(np.concatenate([self.deepgesture_features, self.gene_features, self.snpsift_features], axis=1), dtype=float)
        
        elif not self.run_snpeff and not self.run_dbnsfp:            
            self.dataset = np.array(np.concatenate([self.deepgesture_features, self.gene_features], axis=1), dtype=float)
        else:
            raise Exception('No features were extracted. Please check your input data.')

        snpeff_names = []
        if self.run_snpeff:
            snpeff_names = self.snpeff_names

        snpsift_names = []
        if self.run_dbnsfp:
            snpsift_names += self.snpsift_names

        pickle.dump(
            [
                self.dataset, 
                self.deepgesture_names + self.gene_feature_names  + snpeff_names + snpsift_names, 
                table
            ], 
            open("{}/mutation_features.pk".format(self.output_path), 'wb'), 
            protocol=pickle.HIGHEST_PROTOCOL
        )

        logger.info('[features, feature_names, dataset] stored in: {}/mutation_features.pk'.format(self.output_path))