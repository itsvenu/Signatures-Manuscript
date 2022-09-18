'''
author: @itsvenu_

Input:
    1. Mutational catalogue
    2. COSMIC weights file generated for each denovo profile
    3. Signature reference db
    4. max_iter arg
    5. output folder to save results
'''

import pandas as pd
import signatureanalyzer as sa
import click
import sys, re, os

def format_spectra(input_catalogue, context_type):
    spectra_df = pd.read_csv(input_catalogue, sep = '\t')
    spectra_df = spectra_df.set_index('MutationType')

    if(str(context_type) == '1536'):
        spectra_df = spectra_df.groupby(spectra_df.index.str[1:8]).sum()

    elif(str(context_type) == '83'):
        spectra_df = spectra_df

    else:
        print('Please provide either 1536 or 83 context type')

    return(spectra_df)

## format required COSMIC signatures to list
def cosmic_list(req_cosmic):

    req_sigs = pd.read_csv(req_cosmic, sep = '\t', names = ['sig'])
    req_sigs = req_sigs.sig.tolist()
    return(req_sigs)

def format_refdb(ref_signatures):

    ref_db = pd.read_csv(ref_signatures, sep = '\t')
    ref_db = ref_db.set_index('Type')
    return(ref_db)

## main function
def SupervisedARDNMF(input_catalogue, req_cosmic, ref_signatures, max_iter, context_type, output):

    ## helper functions

    spectra_df = format_spectra(input_catalogue, context_type)
    Wref = format_refdb(ref_signatures)
    req_sigs = cosmic_list(req_cosmic)

    Wref_df = Wref.loc[:, req_sigs]

    ## supervised ARD-NMF
    res_supervised = sa.supervised_bnmf.supervised_ardnmf(spectra_df,
                                                          Wref_df,
                                                          objective='poisson',
                                                          verbose=True,
                                                          max_iter=max_iter)

    ## write. output files
    r_expo = res_supervised['Hraw'].T
    r_expo = r_expo.reset_index(level = 0).rename(columns = {"index": "PID"})
    r_expo_file = output + '/SupervisedExposures.txt'
    r_expo.to_csv(r_expo_file, index=None, sep='\t')

    r_sigs = res_supervised['Wraw']
    r_sigs = r_sigs.reset_index(level = 0).rename(columns = {"index": "MutationType"})
    r_sigs_file = output + '/SupervisedSignatures.txt'
    r_sigs.to_csv(r_sigs_file, index=None, sep='\t')


## command line args
@click.command()
@click.option('--input_catalogue',
              type=click.Path(exists=True),
              help="input mutational catalogue, SA compatible",
              required=True)

@click.option('--req_cosmic',
              help="required COSMIC signatures for supervised ARD-NMF",
              required=True)

@click.option('--ref_signatures',
              type=click.Path(exists=True),
              help="COSMIC reference signatures, either SBS96 or ID83",
              required=True)

@click.option('--max_iter',
              type=click.INT,
              help="maximum number of iteration to run ARD-NMF",
              required=True)

@click.option('--context_type',
              help="context type e.g. 1536, 96, 83",
              required=True)

@click.option('--output',
              help="output path to store results",
              required=True)

##
def all_run(input_catalogue, req_cosmic, ref_signatures, max_iter, context_type, output):

    SupervisedARDNMF(input_catalogue,
                     req_cosmic,
                     ref_signatures,
                     max_iter,
                     context_type,
                     output)

if __name__ == '__main__':
    all_run()

## END ##
