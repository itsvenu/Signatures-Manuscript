'''
SigProfiler - decomposition of denovo signatures to latest COSMIC signatures
'''

import sys, re, os
import pandas as pd
import click
from SigProfilerAssignment import Analyzer as Analyze

def SigProfilerDecomposition(input_catalogue, denovo_signatures, reference_signatures, output):

    # prepare reference signatures
    ref = pd.read_csv(reference_signatures, sep = '\t')
    ref = ref.set_index('Type')

    Analyze.decompose_fit(samples = input_catalogue,
                          output = output,
                          signatures = denovo_signatures,
                          signature_database = ref,
                          verbose = True)


## command-line args
@click.command()
@click.option('--input_catalogue',
              type=click.Path(exists=True),
              help="input mutational catalogue",
              required=True)

@click.option('--denovo_signatures',
              type=click.Path(exists=True),
              help="denovo signatures file",
              required=True)

@click.option('--reference_signatures',
              type=click.Path(exists=True),
              help="COSMIC reference signatures",
              required=True)

@click.option('--output',
              help="output path to store results",
              required=True)

def all_run(input_catalogue, denovo_signatures, reference_signatures, output):
    SigProfilerDecomposition(input_catalogue, denovo_signatures, reference_signatures, output)

if __name__ == '__main__':
    all_run()
