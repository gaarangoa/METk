import pickle
import pandas as pd
import numpy as np 
import os
import subprocess
import logging
import sys

from .parser import variants_parser, variants_parser_large
from .vectors import compute_embeddings
from .parser import format_input_data

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(logging.StreamHandler(stream=sys.stdout))

ROOT_DIR, _ = os.path.split(__file__)

class EmbedDoc():
    def __init__(self, model_name):
        ''' 
        This Class is a wrapper for the C++ implementation of document embeddings on StarSpace to embed a document and retrieve the document embeddings. 
        ...
        
        Attributes
        ----------

        model_name: str
            Path to the pretrained model with StarSpace

        
        Methods
        ----------
        load()
            Loads the pretrained StarSpace model

        embed(input_string='')
            Returns the embeddings from the pretrained model.

        '''
        self.model_name = model_name
        self.load()
    
    def load(self):
        ''' 
        Parameters
        ----------
        '''
        self.model = subprocess.Popen(
            [ '{}/embed_doc_individual.bin'.format(ROOT_DIR), self.model_name],
            shell=False,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            bufsize=0
        )

        for _ in range(3):
            log.info(self.model.stdout.readline().decode("utf-8").strip())

    def embed(self, input_string=''):
        ''' 
        Parameters
        ----------
        input_string : str
            The input string with the tokens to be processed by the model. 

        Returns
        ----------    
        numpy array 
            It returns a numpy array of the dimension of 
            the embeddings defined in the pretrained model as a numpy array.
        '''
        code = self.model.stdin.write(str.encode("{}\n".format(input_string) ))        
        return np.array(self.model.stdout.readline().decode("utf-8").strip().split(), dtype=float)
    
    def close(self, ):
        self.model.kill()

class Embedder():
    def __init__(self, path_to_reference_genomes='', path_to_trained_embeddings='', cores=10, pretrained_embeddings=''):
        
        ''' 
            path_to_reference_genome: Path where the reference genomes are stored
            path_to_pretrained_embeddings: Path where the pretrained embeddings are stored
            pretrained_embeddings: List of pretrained embeddings to run on deepgesture. Default 128 embeddings.
        '''
        
        self.variants_column_names = ",".join(['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Allele_1', 'Tumor_Allele_2', 'Variant_Type', 'Hugo_Symbol'])
        self.filter_mutations = ''
        self.reference_genomes = path_to_reference_genomes
        self.model_path = path_to_trained_embeddings
            
        self.no_embeddings = False
        self.cores = cores
        self.pretrained_embeddings = pretrained_embeddings
        self.debug = False

        
    def call(self, table, reference_genome, identifier='', output_path='', selected_variant_types='INS,DEL,SNP'):
        '''
            table: Input formatted pandas dataframe
            reference_genome: reference genome used in the analysis
            identifier: Column name that will be used as ID for getting the embeddings. If the identifier corresponds to 
                        sample ids. It will calculate embeddings at the patient level. If identifier are unique ids for each 
                        entry, it will calculate embeddings at the mutation level. The identifier can be used to group 
                        entries in any form. 
            output_path: Path where to store the output of deepgesture
            selected_variant_types: A string with the names of the type of variants to process. As default, it uses 
                        INS,DEL,SNP. As the names can change you need to setup the right names. For instance, for ctDNA data from
                        Guardant Health the names are small insertion,small deletion,Single nucleotide substitution. The order has to 
                        be Insertions,deletions,SNVs.
        '''
        
        os.system('mkdir -p {}'.format(output_path))
        
        self.output_path = "./" if not output_path else output_path
        self.variants_reference_genome = reference_genome
        
        self.variants = "{}/deepgesture.tsv".format(output_path)
        self.sample_id_column = identifier
        self.variant_types = selected_variant_types
        
        table.to_csv( "{}".format(self.variants), sep='\t', index=False )
        table.to_vcf( "{}/deepgesture.vcf".format(output_path) )
        
        process(self)
    
def load_embeddings(embeddings_file):
    ''' 
        INPUT:
            embeddings_file: Load embeddings from this pretrained embedder.
            
        RETURNS:
            vectors: Embeddings in a numpy array.
            sample_ids: corresponding ID of the vectors.
    '''
    
    sample_ids = np.array([i.strip() for i in open("{}/samples.txt".format( "/".join(embeddings_file.split('/')[:-1]) ))])
    vectors = np.array([i.strip().split() for i in open("{}".format(embeddings_file))])

    return vectors, sample_ids

def process(args):
    if args.cores > 1:
        args.cores = args.cores - 1
    
    if args.debug == False: 
        logger.setLevel(logging.INFO)

    # make output directory
    args.output_path = args.output_path
    if not os.path.exists( "{}".format(args.output_path) ):
        os.makedirs("{}".format(args.output_path))

    # parse vcf/mutation file
    variants_parser_large(args)

    # obtain embeddings from variants
    if args.no_embeddings == True:
        logger.warning('Embeddings are not being computed, if this is not desired, please remove --no_embeddings')
        return 0
    
    for model, model_name in [i.split(':') for i in args.pretrained_embeddings.split(',')]:
        compute_embeddings(
            input_file='{}/tokens.txt'.format(args.output_path), 
            outdir=args.output_path, 
            model='{}/{}'.format(args.model_path, model),
            model_name=model_name
        )

def compute_only_embeddings(args):
    for model, model_name in [i.split(':') for i in args.pretrained_embeddings.split(',')]:
        compute_embeddings(
            input_file='{}'.format(args.input), 
            outdir=args.output_path, 
            model='{}/{}'.format(args.model_path, model),
            model_name=model_name,
            prefix=args.prefix
        )

