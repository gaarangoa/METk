import pandas as pd
import numpy as np 
import logging 
import sys
from tqdm import tqdm
import json
from pyfaidx import Fasta
import pickle as pk
import os 
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import itertools
import multiprocessing as mp
import math

logger = logging.getLogger()

ROOT_DIR, _ = os.path.split(__file__)

# Variant Parser
def variants_parser():
    pass

def iupac(key):
    codes = {
        ('A','G'): 'R',
        ('C','T'): 'Y',
        ('C','G'): 'S',
        ('A','T'): 'W',
        ('G','T'): 'K',
        ('A','C'): 'M',
        ('C', 'G', 'T'): 'B',
        ('A', 'G', 'T'): 'D',
        ('A', 'C', 'T'): 'H',
        ('A', 'C', 'G'): 'V'
    }

    try:
        return codes[key]
    except:
        return key[1]

def update_sample(samples, _id, key, value):
    try:
        samples[_id][key].append(value)
    except:
        samples[_id].update({key: value})

class ParseVariants(object):
    def __init__(self, args):
        self.args = args
        self.given_column_names = args.variants_column_names.split(',')
        self.bp_offset = 10
        self._set_kmers = [3,4,5]

        # if args.debug:
        logger.setLevel(logging.DEBUG)
        
        logger.info('Processing: {} variants file'.format(args.variants))
        logger.info('Reference genome: {}'.format(args.variants_reference_genome))

        if args.filter_mutations:
            self.variant_classification_column, self.value_to_filter = args.filter_mutations.split(':')

        self.variant_types = args.variant_types.split(',')

        logger.info("Columns used: {}".format(args.variants_column_names) )
        logger.info("Consequence types: {}".format(self.variant_types))

        self.hg = {
            'GRCh37': Fasta('{}/GRCh37.d1.vd1.fa'.format(self.args.reference_genomes)), 
            'GRCh38': Fasta('{}/GRCh38.d1.vd1.fa'.format(self.args.reference_genomes))
        }

    def load_data(self, sample_id=''):
        logger.info('Parsing input file {}'.format(self.args.variants))
        self.data_dict = {}
        for ix,i in tqdm(enumerate(open(self.args.variants, mode='r', encoding='utf8'))):
            if i[0] == '#': continue
            if ix == 0:
                self.column_names = i.strip().split('\t') + ['context', 'iupac_context']
                ft = open('{}/metadata_header.txt'.format(self.args.output_path), 'w')
                ft.write(i)
                ft.close()

                continue

            i = i.strip().split('\t')
            item = {self.column_names[kx]: k for kx,k in enumerate(i)}

            try:
                self.data_dict[item[self.args.sample_id_column]].append(item)
            except:
                self.data_dict[item[self.args.sample_id_column]] = [item]
        
        logger.info("Total entries in the input file: {}".format(len(self.data_dict)))

    def process_core(self, ):
        
        hg = {'GRCh37': Fasta('{}/GRCh37.d1.vd1.fa'.format(self.args.reference_genomes)),  'GRCh38': Fasta('{}/GRCh38.d1.vd1.fa'.format(self.args.reference_genomes))}

        ft = open('{}/tokens.txt'.format(self.args.output_path), 'w', encoding='utf8')
        fm = open('{}/metadata.txt'.format(self.args.output_path, ), 'w', encoding='utf8')
        fs = open('{}/samples.txt'.format(self.args.output_path, ), 'w', encoding='utf8')

        fm.write('{}\n'.format( "\t".join(self.column_names) ))

        filtered_ix = 0
        for ix, sample in enumerate(self.data_dict.values()):
            try:
                sample_value = process_mutations(sample, self.given_column_names, self.args, hg, self.bp_offset, self.variant_types, self._set_kmers, fm)
            except Exception as inst:
                # log.error('Error during execution: {}'.format(inst))
                assert(False)

            if sample_value:
                ft.write("\t".join([" ".join(i) for i in  sample_value]) + '\n')
                fs.write("{}\t{}\n".format(sample[0][self.args.sample_id_column], True))
            else: 
                filtered_ix += 1
                fs.write("{}\t{}\n".format(sample[0][self.args.sample_id_column], False))
        
        logger.info("Processing {} entries in core".format(ix))
        logger.info("Filtered entries {} in core ".format(filtered_ix))

        ft.close()
        fm.close()
        fs.close()

    def chunks(self, p=10):
        it = iter(self.data_dict.values())
        n = math.ceil( len(self.data_dict) / (p+1) )
        for ix,i in enumerate(range(0, len(self.data_dict), n)):
            yield tuple(itertools.islice(it, n)), ix, self.given_column_names, self.args, self.bp_offset, self.variant_types, self._set_kmers
    
    def multiprocessing(self, func, args, workers):
        with ThreadPoolExecutor(workers+1) as ex:
            res = ex.map(func, args)

    def main(self):
        chunks = self.chunks(p=self.args.cores)
        self.multiprocessing(cpu_process, chunks, workers=self.args.cores)

    def merge(self):
        for key in ['tokens', 'metadata', 'samples']:
            
            tokens = " ".join(["{}/{}_{}.txt".format(self.args.output_path,key ,i) for i in range(self.args.cores + 1)])
            if key == 'metadata':
                tokens = "{}/metadata_header.txt ".format(self.args.output_path) + " ".join(["{}/{}_{}.txt".format(self.args.output_path,key ,i) for i in range(self.args.cores + 1)])

            os.system("cat {} > {}/{}.txt".format(tokens, self.args.output_path, key))
            # os.system("rm {}".format(tokens))
        
        logger.info('Merging files: ')
   
def process_mutations(selected_sample, given_column_names, args, hg, bp_offset, variant_types, _set_kmers, fm):
    ''' Process one sample at the time. If the sample is composed of patient data
        it processes all mutations in the patient. If the sample corresponds to 
        individual mutations. It process one mutation at the time. 
    '''
    sample_value = []
    chr_column_name, start_column_name, end_column_name, reference_column_name, variant_1_column_name, variant_2_column_name, variant_type_column_name, gene_column_name = given_column_names
    filtered_entries = 0

    # i are mutations in selected sample
    for ix,i in enumerate(selected_sample):
        gene, variant_type, chrom, start, end, reference, variant_1, variant_2  = i[gene_column_name], i[variant_type_column_name], i[chr_column_name], i[start_column_name], i[end_column_name], i[reference_column_name], i[variant_1_column_name], i[variant_2_column_name]
        start, end = int(start), int(end)

        if reference == variant_1:
            variant = variant_2
        else:
            variant = variant_1

        if args.sample_id_column:
            sample_id = i[args.sample_id_column]
        else:
            sample_id = 'sample_i'

        # logger.debug((gene, variant_type, chrom, start, end, reference, variant))
        chrom = chrom if 'chr' in str(chrom) else 'chr{}'.format(chrom)
        
        if str(gene) == 'nan':
            gene = ''
        
        try:
            seq = hg[args.variants_reference_genome][chrom][(start-bp_offset-1):(end+bp_offset)]
        except:
            logger.debug('Filtering out {} because of no sequence found'.format(i))
            filtered_entries += 1
            # sample_value.append([])
            continue

        # logger.debug(('SeqInfo', seq, len(seq), variant_types, variant_type))

        var_seq = ''
        if variant_type == variant_types[1]: # DEL

            # print(seq, seq[0:bp_offset], seq[-1*bp_offset:])
            var_seq = seq[0:bp_offset].seq + seq[-1*bp_offset:].seq

        elif variant_type == variant_types[0]: # INS
            # This is done for the newest version. 
            seq = hg[args.variants_reference_genome][chrom][(start-bp_offset-1):(start+bp_offset)]
            # print(seq, inseq, seq[0:bp_offset+1], inseq, seq[-1*bp_offset-1:])
            var_seq = seq[0:bp_offset].seq + variant + seq[-1*bp_offset-1:].seq

        elif variant_type == variant_types[2]: # SNV
            inseq = [reference, variant]
            inseq.sort()

            var_seq = seq[0:bp_offset].seq + iupac((*inseq,)) + seq[-1*bp_offset:].seq

        else:
            # sample_value.append([])
            continue

        if not var_seq:
            logger.debug('Filtering out {} No sequence found'.format(i))
            filtered_entries += 1
            # sample_value.append([])
            continue

        kmers = []
        for k in _set_kmers:
            k = int(k)
            kmers += ["".join(list(var_seq)[i:i+k]) for i in range(0, len(var_seq)-k+1, 1)]
        
        kmers = " ".join(kmers)

        sample_value.append(
            [kmers, gene, '{}>{}'.format(reference, variant)]
        )
        
        i['context'] = seq
        i['iupac_context'] = var_seq

        fm.write( "\t".join([str(k) for k in i.values()]) + '\n' )
    
    return sample_value

def variants_parser_large(args):
    parser = ParseVariants(args)
    parser.load_data()
    # parser.main()
    # parser.merge()
    parser.process_core()

# INPUT DATA
def format_input_data(args):
    V=301
    tcga_embeddings = pd.read_csv(
        '{}/{}.wv'.format(args.output_path, 'tcga'),
        names=["V{}".format( str(i) ) for i in range(V)],
        sep='\t', 
    )

    cbio_embeddings = pd.read_csv(
        '{}/{}.wv'.format(args.output_path, 'cbioportal'),
        names=["V{}".format( str(i) ) for i in range(V)],
        sep='\t', 
    )

    del tcga_embeddings['V{}'.format(V-1)]
    del cbio_embeddings['V{}'.format(V-1)]

    cnvs = json.load(open('{}/cnv.maps.json'.format(args.output_path)))
    sample_ids = [i.strip() for i in open('{}/samples.txt'.format(args.output_path))]

    # comments = open('{}/info.txt'.format(args.output_path), 'w')
    # comments.write('sample_id\thas_errors\tmessage\n')

    cnv_maps = []
    sample_data = []
    for p in sample_ids:
        try:
            cnv_maps.append(cnvs[p])
            # comments.write('{}\t0\t-\n'.format(p))
            sample_data.append({'sample_id': p, 'has_errors': 0, 'message': '-'})
        except:
            cnv_maps.append(np.zeros(99))
            # comments.write('{}\t1\tsample does not have CNV data\n')
            sample_data.append({'sample_id': p, 'has_errors': 1, 'message': 'Sample does not have CNV data.'})
            logger.warning('Sample {} does not have CNV data. Adding CNV as all zeros. Please check your input!'.format(p))

    # comments.close()

    cnv_maps = np.array(cnv_maps)
    cnv_maps = pd.DataFrame(cnv_maps)
    cnv_maps = cnv_maps.fillna(0)

    logging.info('TCGA embeddings: {}'.format(tcga_embeddings.shape))
    logging.info('CBIOPORTAL embeddings: {}'.format(cbio_embeddings.shape))
    logging.info('CNV distributions: {}'.format(cnv_maps.shape))

    return tcga_embeddings, cbio_embeddings, cnv_maps, sample_data