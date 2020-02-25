__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""


import datetime
import sys
import os
import pandas as pd
import json

CONTRASTS,=glob_wildcards('input/{contrast}.diffexp.tsv')

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

for contrast in CONTRASTS:
    message("Contrast " + contrast + " will be processed")

rule all:
    input:
        expand(["results/{contrast}/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt", "results/{contrast}/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"], contrast=CONTRASTS,FC=config["FC"], adjp=config["adjp"])


include: "rules/seurat.smk"
