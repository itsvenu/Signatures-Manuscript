# @itsvenu_
# ''' Calculate optimal contribution of mutations to PCAWG SBS/ID signatures'''

import sys, re, os
import click
import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze

def ComputeAttributions(input_catalogue, output, reference_signatures):

    Analyze.cosmic_fit(samples = input_catalogue,
                       output = output,
                       signatures=None,
                       signature_database=reference_signatures,
                       genome_build="GRCh37",
                       verbose=True,
                       collapse_to_SBS96=False)

##
@click.command()
@click.option('--input_catalogue',
              type=click.Path(exists=True),
              help="input mutational catalogue",
              required=True)

@click.option('--output',
              help="output path to store results",
              required=True)

@click.option('--reference_signatures',
              type=click.Path(exists=True),
              help="SBS/ID reference signatures",
              required=True)

def all_run(input_catalogue, output, reference_signatures):
    ComputeAttributions(input_catalogue, output, reference_signatures)


if __name__ == '__main__':
    all_run()
