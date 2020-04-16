# Variant_annotation_Tempus

A small software program to output a CVS table and a VCF file with annotation for each variant in the file. Additional information for each variant is gathered via ExAC. To speed up computation, variants are queried in bulk with ExAC API, and  annotation information is stored in list then append to a dataframe. 


### Dependancies:
Language: python 3

Package used: json, argparse, requests, pandas, numpy, os.path

### Annotated information: 
1. Reads depth at the site of variation; 
2. Number of reads supporting the variant; 
3. Reads supporting the variant versus those supporting reference reads; 
4. The Allele Frequency of variant from EXaC;
5. The consequence of variant from EXaC - the most severe consequence is selected when multiple consequences are available;
6. The related GeneID from EXaC - the related GeneIDs are picked based on the selected variant consequence.

### Usage: 

The input is a VCF file, and the output will be "output.csv" and "output.vcf".

      python variant_annoation_VCF.py -i input.vcf -o output


If an output file is not specified, the output will be written to "input_annotated.csv" and "input_annotated.vcf".
                 
      python variant_annoation_VCF.py -i input.vcf



