# written to more efficiently find the gene
# neighborhood. Depth-first search 
# adding augment-network function
# written 6-7-17 JLW

import itertools, os, sys, time, pickle
import networkx as nx
from collections import defaultdict

def pairwise(iterable):
	a,b = itertools.tee(iterable)
	next(b,None)
	return zip(a,b)

def score_path(pth): # pass a list of nodes in a shortest path
	if '@' not in pth:
		return 1.0
	else: # create pairwise edges, pull scores
		split_path = pth.split('@')
		prws = pairwise(split_path)
		edge_scores = [GENE_GRAPH.get_edge_data(a,b)['weight'] for (a,b) in prws]
		score=1
		for e in edge_scores:
			score*=e
		return score

def extend_path(pth): # should be a list of paths, will return a list of new paths which include the initial path
	if '@' in pth:
		split_path = pth.split('@')
		last_gene = split_path[-1]
		former_path = split_path[-2]

		nbrs = [n for n in nx.all_neighbors(GENE_GRAPH,last_gene)]
		new_paths = [pth+'@'+nbr for nbr in nbrs if nbr!=former_path]

	else: # when a signle gene is passed
		last_gene = pth

		nbrs = [n for n in nx.all_neighbors(GENE_GRAPH,last_gene)]
		new_paths = [pth+'@'+nbr for nbr in nbrs]

	return new_paths

def fast_track(pstack,pth,sterm): # method to fast track some paths that are alternate routes to an approved node
	split_path = pth.split('@')	
	lastn = split_path[-1]

	redun = [p for p in pstack if lastn==p.split('@')[-1]]#remove paths that terminate on an included node
#	if len(redun) >0:
#		print('redundant paths found:')
#		print(redun)
	fast_track = set([p for p in redun if score_path(p)>sterm])
	new_stack = list(set(pstack).difference(set(redun))) #remove redundant
	return (new_stack,fast_track) # return trimmed stack, and only the high-scoring redundant paths

def find_neighborhood(gene,sterm): # sterm = termination path score
	stack = [gene]
	approved = set()
	#stime = time.time()
	while stack:
		#print(sys.getsizeof(stack))
		lastp = stack[-1]
		if score_path(lastp)>sterm:
			approved.add(lastp)
			stack+=extend_path(lastp)
		stack.remove(lastp)	
		(stack,fast) = fast_track(stack,lastp,sterm)
		approved = approved.union(fast) 
        # reformat dictionary before returning                          
	ps_d = defaultdict(float)
	for pth in approved:
		ps_d[pth] = score_path(pth)
	return ps_d
	#print("Time Taken: %.3f sec\n" % (time.time() - stime))
