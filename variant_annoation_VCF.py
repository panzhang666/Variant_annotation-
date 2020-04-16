#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 02:32:32 2020

@author: panzhang
"""

import json
import sys
import requests
import argparse
import os.path
import pandas as pd
import numpy as np


#Add arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="File path to VCF file.")
parser.add_argument("-o", type=str, help="File name for output file.")
args = parser.parse_args()


#Check if input VCF file exists
if not args.i:
	print("No input file specified.")
	exit()

elif not os.path.isfile(args.i):
	print("File specified does not exist.")
	exit()


#Creat names for output annotated VCF file and annotated table 
#Check if output file name exists 
if not args.o:
    inVcf = args.i
    name = inVcf.split(".vcf")[0]
    outputVcf = str(name)+"_annotated.vcf"
    outputTable = str(name)+"_annotated.csv"

elif args.o:
    outputVcf = str(args.o)+".vcf"
    outputTable = str(args.o)+".csv"


"""
Store putative impact for some Variant Consequences to help 
annotate variant with the most deleterious possibility.
Values based on http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf:
0 - MODIFIER, 1 - LOW, 2 - MODERATE, 3 - HIGH
"""

consequenceSeverity = {'3_prime_UTR_variant': 0,
                        '5_prime_UTR_variant': 0,
                        'intron_variant': 0,
                        'non_coding_transcript_exon_variant': 0,
                        'splice_region_variant': 1,
                        'synonymous_variant': 1,
                        'stop_retained_variant': 2,
                        'missense_variant': 2,
                        'initiator_codon_variant': 2,
                        'stop_lost': 3,
                        'stop_gained': 3,
                        'splice_donor_variant': 3,
                        'splice_acceptor_variant': 3}


def get_var_list(vcfFile):
    """ 
    Function to loop through the VCF file to create a list that contains the variants in a format that the
    EXaC API accepts. The format is "[\"chr-position-ref-alt\",\"chr-position-ref-alt\"]" .
    Input: vcf file
    Output: API formatted variants
    """   
    
    with open(vcfFile, 'r') as vcf:
        apiFormattedVariants = []
        for line in vcf:
            # Ignore header lines
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, *rest = line.strip().split('\t')
            apiFormattedVariants.append(r'\"' + "-".join((chrom, pos, ref, alt)) + r'\"')
    return(apiFormattedVariants)
    


def get_exac_info(apiFormattedVariants):
    """
    Function to extract information from EXaC for a list of variants.
    Input: A list of API formatted variants
    Output: An array of information about the variants found in the ExAC database
    """
    
    # The /rest/bulk/variant stores all the variant 
    url = 'http://exac.hms.harvard.edu/rest/bulk/variant' 
    data = r'"[' + ",".join(apiFormattedVariants[0:len(apiFormattedVariants)]) + r']"'
    # eval() is required so that the variants are formatted in the way the API expects
    r = requests.post(url, eval(data))
    exacInformation = json.loads(r.text)
    return(exacInformation)


def extract_exac_annot(exacInformation, consequenceSeverity):
    """
    Function to extract helpful annoation from EXaC information.
    Input: 1. A dictionary stores information of variants from EXaC for variants
           2. A dictionary stores putative impact for some Sequence Ontology terms 
    Output: Three dictionaries store "Allele frequency", "Variant Consequence", 
            "Related GeneID" information for variants separately
    """
    
    variantConsequences = {}
    exacAlleleFreq = {}
    variantGeneID = {}
    for variant in exacInformation:
        # Use 'chr-position' as the key for the exac_allel_freq and variantConsequences dictionaries
        position = "-".join(variant.split('-')[0:2])
        # Get the allele frequencies from EXaC
        if 'allele_freq' in exacInformation[variant]['variant']:
            allele_freq = '{:.3f}'.format(exacInformation[variant]['variant']['allele_freq'])
        else:
            allele_freq = 'NA'

        exacAlleleFreq[position] = allele_freq

        # Obtain the variant consequence from EXaC
        if exacInformation[variant]['consequence'] is None or exacInformation[variant]['consequence'] == {}:
            variantConsequences[position] = 'NA'
            variantGeneID[position] = 'NA'
        else:
            # Loop through consequences for the variants that have more than one consequence
            if len(exacInformation[variant]['consequence']) > 1:
                highestSevereVal = -1
                mostSevereConsequence = None

                for consequence in exacInformation[variant]['consequence']:
                    if consequenceSeverity[consequence] == highestSevereVal:
                        # If there are multiple consequences of the same level keep the one
                        # that sorts first alphabetically
                        if consequence < mostSevereConsequence:
                            mostSevereConsequence = consequence

                    if consequenceSeverity[consequence] > highestSevereVal:
                        highestSevereVal = consequenceSeverity[consequence]
                        mostSevereConsequence = consequence

                variantConsequences[position] = mostSevereConsequence
                
                #Get all the related geneID corresponding to mostSevereConsequence
                #If there is more than one related geneID, join all of them by "," 
                if exacInformation[variant]['consequence'][mostSevereConsequence] is None or exacInformation[variant]['consequence'][mostSevereConsequence] == {}:                    
                    variantGeneID[position] = 'NA'
                else:
                    variantGeneID[position] = ",".join(exacInformation[variant]['consequence'][mostSevereConsequence])                   
            else:
                # Get the consequence for variants that only have one
                variantConsequences[position] = list(exacInformation[variant]['consequence'].keys())[0]
                variantGeneID[position] = ",".join(exacInformation[variant]['consequence'][variantConsequences[position]])
                continue
            
    return(exacAlleleFreq, variantConsequences, variantGeneID)

    
def get_basic_info(vcfLine):
    """
    Function to retrieve relevant annotations from each variant hit
    Input: Line in VCF file
    Output: Tab separated annotations
    """    
    
    chrom, pos, ID, ref, alt, qual, Filter, info, Format, *genotypes = vcfLine.strip().split('\t')
    
    #Create dictionary of info field, value pairs from the INFO section
    separateInfo = info.split(";")
    infoDict = dict([x.split('=') for x in separateInfo])

    #Type of variation 
    typeVar = infoDict['TYPE']

    #Depth of sequence coverage at the site of variation
    totalDepth = infoDict["DP"]

    #Number of reads supporting the variant
    altDepth = infoDict["AO"]
    
    refDepth = infoDict["RO"]
    #Percentage of reads supporting the variant 
    try: 
        if ',' in altDepth:
            varPercent = []
            for depth in altDepth.split(','):
                varPercent.append('{:.3f}'.format((float(depth)/float(refDepth))))
            varPercent = ','.join(varPercent)
        else:
            varPercent = '{:.3f}'.format(float(altDepth)/float(refDepth))

    except ZeroDivisionError: #Instance when there are zero reads supporting the reference
        varPercent = "NA"
    
    extractLine = "\t".join([chrom, pos, ID, typeVar, ref, alt, totalDepth, altDepth, varPercent])
    return (extractLine, totalDepth, altDepth, varPercent )


"""
Function to make a flat list out of list of lists.
"""
flatten = lambda l: [item for sublist in l for item in sublist]


"""
Write the annoation to two files: an annotated VCF file and 
a csv file includes basic variant information and annotation.
"""

vcfFile = args.i  
outputVcf = open(outputVcf, 'w')

#Get a list of EXaC API formatted variants from the input VCF file
apiFormatVariantsList = get_var_list(vcfFile)

#Get an array of information about the input variants found in the ExAC database.
exacInformation = get_exac_info(apiFormatVariantsList)

# Get three dictionaries store "Allele frequency", "Variant Consequence", 
# and "Related GeneID" information for variants separately
exacAlleleFreq, variantConsequences, variantGeneID = extract_exac_annot(exacInformation, consequenceSeverity)

#Create dataframe for variants annotation table
output_table = pd.DataFrame()
output_list =[]


#Process VCF file
with open(vcfFile, 'r') as vcf:
        apiFormattedVariants = []
        for line in vcf:
            # Ignore header lines
            if line.startswith("#"):
                print(line, end='', file=outputVcf)
                #Add description for annotation that will be added
                if line.startswith("##INFO=<ID=END"):
                    print('##INFO=<ID=ANNOT,Description="Variant annotations. Total read depth | Num Reads Supporting Variant | Variant reads vs Ref reads | EXaC variant frequency | Variant Consequence | Related GeneID ">', file=outputVcf)
                continue
            
            chrom, pos, ID, ref, alt, qual, Filter, info, Format, *genotypes = line.strip().split('\t')
            info = info.split(';')
            
            #Get the basic annoation information from input VCF file
            vcf_info, totalDepth, altDepth, varPercent = get_basic_info(line)
            position = chrom + "-" + pos
            
            #Combine the annotation from VCF file and EXaC for csv and VCF output 
            total_annotation = flatten([vcf_info.strip().split("\t"), [exacAlleleFreq[position],variantConsequences[position], variantGeneID[position]]])
            output_list.append(total_annotation)
            
            new_info = ';'.join(info) + ';ANNOT=' + ('|').join((totalDepth, altDepth, varPercent, exacAlleleFreq[position], variantConsequences[position], variantGeneID[position]))            
            print(chrom, pos, ID, ref, alt, qual, Filter, new_info, Format, '\t'.join(genotypes), sep='\t', file=outputVcf)

# Store annoation information in the dataframe and then write to a CSV file
output_table = output_table.append(output_list)
output_table.columns= ["Chrom","Pos","ID", "Type", "Ref", "Alt", "Depth","Num Reads Supporting Variant","Variant reads vs Ref reads", "EXaC variant frequency", "Variant Consequence", "Related GeneID"]           
output_table.to_csv(outputTable,index=False)
outputVcf.close()
           
"""
Thank you Tempus bioinformatics team for taking the time to review my code!
-Pan Zhang
"""
