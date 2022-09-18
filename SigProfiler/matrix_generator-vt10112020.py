# @itsvenu_
# generate mutations matrices given a path of vcf dir

import sys, re, os
import click
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# matrices = matGen.SigProfilerMatrixGeneratorFunc("test", "GRCh37", "/Users/ebergstr/Desktop/test",plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)

def genrate_matrix(project_name, input_vcf, plot, seq_info, exome):

    matrices = matGen.SigProfilerMatrixGeneratorFunc(project_name, "GRCh37", input_vcf, plot=plot, exome=exome, seqInfo=seq_info)

##
@click.command()
@click.option("--project_name",
              help="name for the output matrix files",
              required=True)

@click.option('--input_vcf',
              type=click.Path(exists=True),
              help="input path of vcf files",
              required=True)

@click.option('--plot',
              is_flag=True,
              default=False,
              help="add --plot if the per-sample profile plots are needed")

@click.option('--seq_info',
              is_flag=True,
              default=False,
              help="add --seq_info to write each mutations context into a file")

@click.option('--exome',
              is_flag=True,
              default=False,
              help="add --exome if the mutational catalogues are from exome sequencing")

##
def all_run(project_name, input_vcf, plot, seq_info, exome):
    genrate_matrix(project_name, input_vcf, plot, seq_info, exome)


if __name__ == '__main__':
    all_run()
