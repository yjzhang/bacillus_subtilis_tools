from collections import defaultdict
import csv
import json
import sqlite3
import pandas as pd
from scipy import sparse
from scipy.io import mmread

# 1. load data
# data source: http://subtiwiki.uni-goettingen.de/v3/exports

# 1.0. load gene_names.csv
# locus, name, description, References.Reviews, References.Research papers
gene_names_table = pd.read_csv('gene_names_description', quotechar='"', skipinitialspace=True, quoting=csv.QUOTE_ALL,engine='python')
gene_names = gene_names_table.name

# 1.1. load geneCategories.csv
# category1, category2, cateogry3, category4, category5
# some of these are na
categories_table = pd.read_csv('geneCategories.csv')
categories_isna = pd.isna(categories_table)

# 1.2. load regulations.csv
# the "term" should be regulon + mode
regulations_table = pd.read_csv('regulations.csv')

# 1.3. load interactions.csv
interactions_table = pd.read_csv('interactions.csv')

# 1.4. load operons.csv
# operons_table.genes is a '-' separated list of gene names
operons_table = pd.read_csv('operons.csv')

# create sqlite3 db
conn = sqlite3.connect('subtiwiki/data/subtiwiki.db')
c = conn.cursor()

# create a table representing a term - gene mapping.
#try:
# term_gene:
# term: one of the gene categories or regulons
# gene: gene name
# category: one of "categories", "regulations", or "interactions"
c.execute('CREATE TABLE term_gene(term text, gene text, type text, level int)')
# iterate over nonzero entries of corpus
for i, row in categories_table.iterrows():
    gene = row.gene
    if not pd.isna(row.category1):
        c.execute('INSERT INTO term_gene VALUES (?, ?, ?, ?)', (row.category1, gene, 'category', 1))
    if not pd.isna(row.category2):
        c.execute('INSERT INTO term_gene VALUES (?, ?, ?, ?)', (row.category2, gene, 'category', 2))
    if not pd.isna(row.category3):
        c.execute('INSERT INTO term_gene VALUES (?, ?, ?, ?)', (row.category3, gene, 'category', 3))
    if not pd.isna(row.category4):
        c.execute('INSERT INTO term_gene VALUES (?, ?, ?, ?)', (row.category4, gene, 'category', 4))
    if not pd.isna(row.category5):
        c.execute('INSERT INTO term_gene VALUES (?, ?, ?, ?)', (row.category5, gene, 'category', 5))
#except Exception as e:
#    print('error in gene categories:', str(e))

# load regulations
#try:
for i, row in regulations_table.iterrows():
    if pd.isna(row['mode']):
        term = row.regulon
    else:
        term = row.regulon + ' ' + row['mode']
    gene = row.gene
    c.execute('INSERT INTO term_gene VALUES (?, ?, ?, ?)', (term, gene, 'regulations', 0))
#except Exception as e:
#    print('error in regulation:', str(e))

# load operons
for i, row in operons_table.iterrows():
    term = row.operon
    genes = row.genes.split('-')
    for gene in genes:
        c.execute('INSERT INTO term_gene VALUES (?, ?, ?, ?)', (term, gene, 'operons', 0))


c.execute('CREATE INDEX term_gene_index on term_gene(term, gene, type, level)')

# create a gene info table
c.execute('CREATE TABLE gene_info(gene text, locus text, description text, reviews text, research_papers text)')
for i, row in gene_names_table.iterrows():
    gene = row['name']
    print(gene)
    locus = row.locus
    desc = row.description
    reviews = row['References.Reviews']
    if pd.isna(reviews):
        reviews = ''
    else:
        reviews = reviews.strip('</pubmed>')
    papers = row['References.Research papers']
    if pd.isna(papers):
        papers = ''
    else:
        reviews = reviews.strip('</pubmed>')
    c.execute('INSERT INTO gene_info VALUES (?, ?, ?, ?, ?)', (gene, locus, desc, reviews, papers))

c.execute('CREATE INDEX gene_info_index ON gene_info(gene)')

conn.commit()
conn.close()
