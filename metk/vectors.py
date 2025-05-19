import subprocess
import logging
import os 

ROOT_DIR, _ = os.path.split(__file__)
ROOT_DIR = ROOT_DIR.replace('/vectors', '/bin')
logger = logging.getLogger()

def compute_embeddings(input_file='', outdir='', model='', model_name='model', prefix=''):

       if prefix:
              prefix  = prefix + '.'

       wv_file = '{}/{}{}.wv'.format(outdir, prefix, model_name)
       logger.info('Processing {}'.format(input_file))
       logger.info('Text model: {}'.format(model))
       logger.info('Storing embeddings to {}'.format(wv_file))

       out = subprocess.Popen(['{}/embed_doc'.format(ROOT_DIR), model, input_file, wv_file], 
              stdout=subprocess.PIPE, 
              stderr=subprocess.STDOUT)
       stdout, stderr = out.communicate()

       logger.info('STDOUT: {}'.format(stdout))
       logger.info('STDERR: {}'.format(stderr))