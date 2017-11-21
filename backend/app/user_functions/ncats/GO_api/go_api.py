import networkx as nx
import numpy as np
import pandas as pd
import goenrich
import os
from copy import deepcopy
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt

import pickle

import goenrich.export

GO_path = os.path.abspath(os.path.dirname(__file__))

DB_path = os.path.join(GO_path,"GO_DB")

"""
GO_api

Runs the Q2 pipeline with the given drug and disease

Input:
QDrug - String of the drug common name
QDisease - String of the disease common name
gen_image - Whether you would like the script to render a network graphic

Output:
results - Pandas dataframe of GO Enrichment for the drug-disease combo with the following columns:

"""
class GO_api:


    def __init__(self, DB_path):

        # Code to regenerate pkl files from original ontology data
        # self.O = goenrich.obo.ontology(os.path.join(DB_path, 'go.obo'))
        # self.gene2go = goenrich.read.gene2go(os.path.join(DB_path, 'gene2go.gz'))
        # self.values = {k: set(v) for k, v in self.gene2go.groupby('GO_ID')['GeneID']}

        pkl_file = open(os.path.join(DB_path, 'GO.pkl'), 'rb')
        self.O = pickle.load(pkl_file)
        pkl_file.close()


        pkl_file = open(os.path.join(DB_path, 'gene2go.pkl'), 'rb')
        self.gene2go = pickle.load(pkl_file)
        pkl_file.close()

        pkl_file = open(os.path.join(DB_path, 'values.pkl'), 'rb')
        self.values = pickle.load(pkl_file)
        pkl_file.close()


    def calc_GO_enrichment(self, query, plotFile, target_list=[]):
        tempO = deepcopy(self.O)
        tempVals = deepcopy(self.values)
        background_attribute = 'gene2go'
        # propagate the background through the ontology

        goenrich.enrich.propagate(tempO, tempVals, background_attribute)

        df = self.analyze(tempO, query, background_attribute, target_list=target_list, gvfile=plotFile + '.dot')
        return(df)

    def get_GO_terms(self, gene):
        terms = set()
        for g in gene:
            temp = self.gene2go.loc[self.gene2go['GeneID'] == g]
            [terms.add(x) for x in temp['GO_ID'].values.tolist()]
        return(list(terms))



    def analyze(self, O, query, background_attribute, target_list, **kwargs):
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
                temp = sig.tolist() + target_list
            else:
                raise NotImplementedError(show)
            # All
            # G = goenrich.enrich.induced_subgraph(O, temp)

            # Drug targets
            # G = goenrich.enrich.induced_subgraph(O, target_list)

            # Drug and disease tight overlap
            # overlap = set(sig.tolist()).intersection(set(target_list))
            # G = goenrich.enrich.induced_subgraph(O, overlap)

            # all_sig overlap
            sig = df.query('q<0.05')['term']
            overlap = set(sig.tolist()).intersection(set(target_list))
            G = goenrich.enrich.induced_subgraph(O, overlap)

            for term, node, q, x, n, rej in zip(terms, nodes, qs, xs, ns, rejs):
                if term in G:
                    G.node[term].update({'name' : node['name'], 'x' : x,
                        'q' : q, 'n' : n, 'significant' : rej, 'drug_target': False})
            self.to_graphviz(G.reverse(copy=False), target_list, **options)
        return df


    def to_graphviz(self, G, target_list, gvfile, graph_label='', **kwargs):
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
                    attr['fillcolor'] = 'yellow'
                    attr['style'] = 'filled'
                    attr['label'] = "{name}\\n{x} / {n} genes\\nq = {q:E}".format(name=node['name'],
                                                                                  q=node['q'], x=node['x'], n=node['n'])
                elif node['isTarget']:
                    attr['fillcolor'] = 'green'
                    attr['style'] = 'filled'
                    attr['label'] = "{name}\\n(drug target)".format(name=node['name'])
                else:
                    attr['fillcolor'] = 'red'
                    attr['style'] = 'filled'
                    attr['label'] = "{name}\\n{x} / {n} genes\\nq = {q:E}".format(name=node['name'],
                                                                              q=node['q'], x=node['x'], n=node['n'])
            else:
                if node['isTarget']:
                    attr['fillcolor'] = 'green'
                    attr['style'] = 'filled'
                    attr['label'] = "{name}\\n(drug target)".format(name=node['name'])
                else:
                    attr['color'] = 'black'
                    attr['fillcolor'] = 'black'
                    attr['height'] = 0.25
                    attr['fixedsize'] = 'true'
                    attr['style'] = 'filled'
                    attr['shape'] = 'circle'
                    attr['label'] = ""

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
    GO_test = GO_api(DB_path)

    disease_gene_list = {7555, 22796, 7182, 3990, 80150, 1071, 4023, 185, 6462, 319, 27202, 5444, 1482, 336, 6352, 338,
                     3156, 344, 348, 7138, 4323, 4837, 3949, 3952, 2169}
    GO_test.calc_GO_enrichment(disease_gene_list, "test", target_list=["GO:0048521", "GO:0005575"])

