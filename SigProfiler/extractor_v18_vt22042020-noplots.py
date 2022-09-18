# command line wrapper for sigprofiler extractor
# input - mutational catalogue in text format, tab separated
# output - path to store results
# min_sig - minimum number of signatures to extrac
# max_sig - maximum number of signatures to extract
# iteration - 1000 (default)
# seq_type - WGS/WES
# mut_type - SNV/ID


import sys, re, os
import pandas as pd
import click
from SigProfilerExtractor import sigpro as sig

## main function
def extract_signatures(input_catalogue, output, min_sig, max_sig, max_iter, cat_channel, exome):

    """ command line wrapper for sigprofiler extractor"""

    cat_channel_format = str(cat_channel).split()

    if not os.path.exists(output):
        os.makedirs(output)

    ## print some useful info
    if exome:
        seq_type = "WES"

    else:
        seq_type = "WGS"

    print("\nExtracting signatures - " + " " + seq_type + "; " + str(cat_channel) + " channel signatures!\n")

    ## extraction
    sig.sigProfilerExtractor(input_type="text",
                             output=output,
                             input_data=input_catalogue,
                             reference_genome='GRCh37',
			     opportunity_genome='GRCh37',
                             minimum_signatures=min_sig,
                             maximum_signatures=max_sig,
                             nmf_replicates=max_iter,
                             context_type=cat_channel_format,
			     make_decomposition_plots=True,
                             exome=exome)

## command line options

@click.command()
@click.option('--input_catalogue',
              type=click.Path(exists=True),
              help="input mutational catalogue",
              required=True)

@click.option('--output',
              help="output path to store results",
              required=True)

@click.option('--min_sig',
              type=click.INT,
              help="minimum number of signatures to extract",
              default=2,
              required=True)

@click.option('--max_sig',
              type=click.INT,
              help="maximum number of signatures to extract",
              required=True)

@click.option('--max_iter',
              type=click.INT,
              help="maximum number of iterations to perform",
              default=1000,
              required=True)

@click.option('--cat_channel',
              help="one of 96, 1536, DINUC, ID - depending on the mutational catalogue",
              required=True)

@click.option('--exome',
              is_flag=True,
              default=False,
              help="add --exome if the mutational catalogues are from exome sequencing",
              required=True)

def all_run(input_catalogue, output, min_sig, max_sig, max_iter, cat_channel, exome):
    extract_signatures(input_catalogue, output, min_sig, max_sig, max_iter, cat_channel, exome)

if __name__ == '__main__':
    all_run()

####
