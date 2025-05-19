import subprocess
import logging
import os 

ROOT_DIR, _ = os.path.split(__file__)
logger = logging.getLogger()

def compute_embeddings(input_file='', outdir='', model='', model_name='model', prefix=''):
        '''
            INPUT:
            input_file: Input file with tokens.
            outdir: Output directory to store embeddings.
            model: Pretrained model to use for embeddings.
            model_name: Name of the model.
            prefix: Prefix for the output file name.
        '''
        if prefix:
                prefix  = prefix + '.'

        try:
            os.system('chmod +x {}/embed_doc'.format(ROOT_DIR))
        except:
            pass

        wv_file = '{}/{}embeddings.wv'.format(outdir, prefix)
        logger.info('Processing {}'.format(input_file))
        logger.info('Text model: {}'.format(model))
        logger.info('Storing embeddings to {}'.format(wv_file))

        out = subprocess.Popen(['{}/embed_doc'.format(ROOT_DIR), model, input_file, wv_file], 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT)
        stdout, stderr = out.communicate()

        logger.info('STDOUT: {}'.format(stdout))
        logger.info('STDERR: {}'.format(stderr))