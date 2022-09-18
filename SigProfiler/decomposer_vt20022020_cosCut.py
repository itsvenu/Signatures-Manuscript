
# signature - probs
# activities - exposures
# catalogue
# output_path, create dir

import sys, re, os
import pandas as pd
import click
from SigProfilerExtractor import decomposition as decomp

def decompose_denovo(signatures, activitis, mut_catalogue, mut_type, sig_db, cosine_sim, output):

    '''Command line wrapper for sigprofiler extractor decomposition'''

    cat_channel = str(mut_type)

    # create output directory
    if not os.path.exists(output):
        os.makedirs(output)

    decomp.decompose(signatures=signatures,
                     activities=activitis,
                     samples=mut_catalogue,
                     output=output,
                     genome_build="GRCh37",
                     signature_database=sig_db,
                     make_decomposition_plots=True,
                     newsignature_threshold=cosine_sim,
                     verbose=True)

## command line options

@click.command()
@click.option('--signatures',
              type=click.Path(exists=True),
              help="de novo signature probabilities",
              required=True)

@click.option('--activitis',
              type=click.Path(exists=True),
              help="de novo signature activities or exposures",
              required=True)

@click.option('--mut_catalogue',
              type=click.Path(exists=True),
              help="mutational catalogue used for extraction",
              required=True)

@click.option('--mut_type',
              type=click.INT,
              help="mutational context, either 96, 78 or 83",
              required=True)

@click.option('--sig_db',
              type=click.Path(exists=True),
              help="custom signature database",
              required=True)

@click.option('--cosine_sim',
              type=click.FLOAT,
              help="cosine similarity cut off",
              required=True)


@click.option('--output',
              help="output path to store results",
              required=True)

def all_run(signatures, activitis, mut_catalogue, mut_type, sig_db, cosine_sim, output):
    decompose_denovo(signatures, activitis, mut_catalogue, mut_type, sig_db, cosine_sim, output)

if __name__ == '__main__':
    all_run()
