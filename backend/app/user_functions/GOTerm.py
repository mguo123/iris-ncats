from iris import state_types as t
from iris import IrisCommand
from iris import iris_objects
from iris import state_machine as sm
from iris import util as util

import pandas as pd


class GotermEnrichment(IrisCommand):
    title = "compute GO term enrichment on {gene_list}"
    examples = [ "compute GO term enrichment on a list of gene symbols" , "describe the function of a list of genes"]
    help_text = [
        "Gene Ontology (GO) Consortium contains functional and structure annotations for genes. Given a set of genes that are up-regulated under certain conditions, an enrichment analysis will find which GO terms are over-represented using annotations for a list of genes. \nCurrenlty, Iris takes in Gene Symbols (e.g. CTSK) and corrects p-values with FDR"
    ]
    argument_types = {
        "gene_list": t.List("What is the list of gene symbols (e.g. CTSK) you want to analyze?"),
        "n_top": t.Int("How many top results you want to see?")
    }
    argument_help = {
        "gene_list": "A list of gene symbols of interest",
        "n_top": "Number of the top enrichment GO terms you want to see",
    }
    def command(self, gene_list, n_top):
        from gprofiler import GProfiler
        import numpy as np
        gp = GProfiler("")
        r0 = gp.gprofile(gene_list,correction_method=GProfiler.THR_FDR,ordered=True)
        r0 = np.array(r0)
        r0 = r0[r0[:,9]=='MF']
        name_out = r0[0:n_top,-3]
        p_out = r0[0:n_top,2]
        data = np.array([x for x in zip(name_out, p_out)])
        dataframe = pd.DataFrame(data=data, columns=["Name", "p-values"])

        return dataframe

    def explanation(self, results):

        iris_df = iris_objects.IrisDataframe(data=results)

        return [ "These are the enriched molecular function GO terms with corrected p-values", iris_df]

gotermEnrichment = GotermEnrichment()
