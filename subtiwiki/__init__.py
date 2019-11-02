import os
import sqlite3

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

PATH = os.path.dirname(__file__)
DB_DIR = os.path.join(PATH, 'data', 'subtiwiki.db')
MODES = ['all', 'category', 'regulations', 'operons']

@lru_cache(maxsize=None)
def get_all_genes(db_dir=DB_DIR):
    """
    Returns a list of all unique gene symbols.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    C.execute('SELECT DISTINCT gene FROM gene_info')
    results = C.fetchall()
    conn.close()
    return [x[0] for x in results]

def get_gene_info(genes, db_dir=DB_DIR, return_header=True):
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    results = []
    for gene in genes:
        C.execute('SELECT * FROM gene_info WHERE gene=?', (gene,))
        r = C.fetchall()
        print(gene, r)
        results += r
    conn.close()
    if return_header:
        results = [('Gene', 'Locus', 'Description', 'Review PMIDs', 'Paper PMIDs')] + results
    return results

@lru_cache(maxsize=None)
def get_all_term_id_names(db_dir=DB_DIR, mode='all'):
    """
    Returns a list of all unique tuples (type, term), where
    type is one of 'category', 'regulations', or 'operons'
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    if mode == 'all':
        C.execute('SELECT DISTINCT type, term FROM term_gene')
    else:
        C.execute('SELECT DISTINCT type, term FROM term_gene WHERE type=?', (mode,))
    results = C.fetchall()
    conn.close()
    return results

@lru_cache(maxsize=None)
def get_term_genes(term, threshold=3, db_dir=DB_DIR):
    """
    Given a term, this returns a list of all genes associated with that term.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    C.execute('SELECT gene FROM term_gene WHERE term=?', (term, ))
    results = C.fetchall()
    results = [x[0] for x in results]
    conn.close()
    return results


def hypergeometric_test(genes, return_header=False, mode='all', db_dir=DB_DIR,):
    """
    Uses a hypergeometric test to identify the most relevant terms.

    Args:
        genes (list): upper-case gene names
        return_header (bool): if True, returns a tuple as the first element.
        mode (str): one of 'all', 'category', 'regulations', 'operons'
        db_dir (str): either DB_DIR, ANATOMY_DIR

    Returns:
        list of 5-tuples: Term Type, Term, p-value, overlapping genes, in order
        of ascending p-value.
    """
    from scipy import stats
    # TODO: make sure genes are correct case - "reverse caps"?
    genes = set(genes)
    all_terms = get_all_term_id_names(db_dir=db_dir, mode=mode)
    all_genes = [x for x in get_all_genes(db_dir=db_dir)]
    cell_p_vals = {}
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    for term_type, term in all_terms:
        genes_term = set(x for x in get_term_genes(term, db_dir=db_dir))
        overlapping_genes = genes.intersection(genes_term)
        if len(overlapping_genes) == 0:
            continue
        k = len(overlapping_genes)
        pv = stats.hypergeom.cdf(k - 1, len(all_genes), len(genes_term), len(genes))
        overlapping_genes = list(overlapping_genes)
        cell_p_vals[term] = (term_type, term, 1 - pv, overlapping_genes)
    cell_p_vals = list(cell_p_vals.items())
    cell_p_vals.sort(key=lambda x: x[1][2])
    # merge items
    cell_p_vals = [x[1] for x in cell_p_vals]
    if return_header:
        header = ['Term Type', 'Term', 'P-value', 'Overlapping Genes']
        cell_p_vals = [header] + cell_p_vals
    return cell_p_vals
