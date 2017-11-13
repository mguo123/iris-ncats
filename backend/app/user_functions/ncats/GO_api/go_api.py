import networkx as nx
import numpy as np
import pandas as pd
import goenrich
import os

from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection

import goenrich.export

DB_path = "./ncats/GO_api/GO_DB"

"""
Q2_query

Runs the Q2 pipeline with the given drug and disease

Input:
QDrug - String of the drug common name
QDisease - String of the disease common name
gen_image - Whether you would like the script to render a network graphic

Output:
results - Pandas dataframe of GO Enrichment for the drug-disease combo with the following columns:
'M'
N' 
'n' 
'name' 
'namespace'
'p'
'q'
'rejected'
'term'
'x'

"""
def calc_GO_enrichment(query, plotFile, target_list=[]):
    O = goenrich.obo.ontology(os.path.join(DB_path, 'go.obo'))

    gene2go = goenrich.read.gene2go(os.path.join(DB_path,'gene2go.gz'))

    values = {k: set(v) for k, v in gene2go.groupby('GO_ID')['GeneID']}

    # propagate the background through the ontology
    background_attribute = 'gene2go'
    goenrich.enrich.propagate(O, values, background_attribute)

    df = analyze(O, query, background_attribute, target_list=target_list, gvfile=plotFile + '.dot')
    return(df)

def get_GO_terms(gene):
    terms = set()
    for g in gene:
        temp = gene2go.loc[gene2go['GeneID'] == g]
        [terms.add(x) for x in temp['GO_ID'].values.tolist()]
    return(list(terms))



def analyze(O, query, background_attribute, **kwargs):
    """ run enrichment analysis for query

    :param O: Ontology graph after backgroud was set
    :param query: array like of ids
    :returns: pandas.DataFrame with results
    """
    options = {
            'show' : 'top20'
    }
    options.update(kwargs)
    _query = set(query)
    terms, nodes = zip(*O.nodes(data=True))
    M = len({x for n in nodes for x in n[background_attribute]}) # all ids used
    N = len(_query)
    ps, xs, ns = goenrich.enrich.calculate_pvalues(nodes, _query, background_attribute,
            M, **options)
    qs, rejs = goenrich.enrich.multiple_testing_correction(ps, **options)
    try:
        df = goenrich.export.to_frame(nodes, term=terms, q=qs, rejected=rejs,
            p=ps, x=xs, n=ns, M=M, N=N)
    except:
        return(None)
    if 'gvfile' in options:
        show = options['show']
        if show.startswith('top'):
            top = int(show.replace('top', ''))
            sig = df.sort_values('q').head(top)['term']
        else:
            raise NotImplementedError(show)
        G = goenrich.enrich.induced_subgraph(O, sig)
        for term, node, q, x, n, rej in zip(terms, nodes, qs, xs, ns, rejs):
            if term in G:
                G.node[term].update({'name' : node['name'], 'x' : x,
                    'q' : q, 'n' : n, 'significant' : rej, 'drug_target': False})
        to_graphviz(G.reverse(copy=False), **options)
    return df


def to_graphviz(G, gvfile, target_list, graph_label='', **kwargs):
    """ export graph of signifcant findings to dot file.
    A png can be generated from the commandline using graphviz

    :param G: the graph to be exported
    :param gvfile: file or filepath
    :param graph_label: For empty label pass graph_label=''.
    """
    for n in G:
        node = G.node[n]
        if n in target_list:
            node['isTarget'] = True
        else:
            node['isTarget'] = False
        attr = {}
        attr['shape'] = 'record'
        if not np.isnan(node.get('q', float('NaN'))):

            if node['isTarget'] and node['significant']:
                attr['color'] = 'yellow'
            elif node['isTarget']:
                attr['color'] = 'green'
            elif node['significant']:
                attr['color'] = 'red'
            else:
                attr['color'] = 'black'
            attr['label'] = "{name}\\n{x} / {n} genes\\nq = {q:E}".format(name=node['name'],
                                                                          q=node['q'], x=node['x'], n=node['n'])
        else:
            if node['isTarget']:
                attr['color'] = 'green'
            else:
                attr['color'] = 'black'
            attr['label'] = """{name}""".format(name=node['name'])
        G.node[n].clear()
        G.node[n].update(attr)

    A = nx.drawing.nx_agraph.to_agraph(G)
    A.graph_attr['label'] = graph_label
    A.graph_attr['labelloc'] = 't'

    if hasattr(gvfile, 'write'):
        A.write(gvfile)
    else:
        with open(gvfile, 'w+') as f:
            A.write(f)


# Paths for unit testing
if __name__ == "__main__":
    O = goenrich.obo.ontology('db/go.obo')

    gene2go = goenrich.read.gene2go('db/gene2go.gz')

    values = {k: set(v) for k, v in gene2go.groupby('GO_ID')['GeneID']}

    # propagate the background through the ontology
    background_attribute = 'gene2go'
    goenrich.enrich.propagate(O, values, background_attribute)

    disease_gene_list = {7555, 22796, 7182, 3990, 80150, 1071, 4023, 185, 6462, 319, 27202, 5444, 1482, 336, 6352, 338,
                     3156, 344, 348, 7138, 4323, 4837, 3949, 3952, 2169}
    calc_GO_enrichment(disease_gene_list, "test", target_list=["GO:0048521", "GO:0005575"])

# Paths for package style loading
else:
    print("other path")
    O = goenrich.obo.ontology(os.path.join(DB_path, 'go.obo'))

    gene2go = goenrich.read.gene2go(os.path.join(DB_path,'gene2go.gz'))

    values = {k: set(v) for k, v in gene2go.groupby('GO_ID')['GeneID']}


