# Variant_annotation_Tempus

A small software program to output a table annotating each variant in the file. Additional information for each variant is gathered via ExAC. To speed up computation time, Variants are queried in bulk with ExAC API and  . Only the most deleterious possibility is selected from each variant entry.

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

Language: Perl 5, version 22
Modules used: JSON, LWP, HTTP, Getopt


