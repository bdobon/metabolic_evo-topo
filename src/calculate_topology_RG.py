#!/usr/bin/env python

'''

    Calculate topological measures of a DIRECTED graph. 
    
    Mainly from:
    http://networkx.readthedocs.org/en/stable/reference/algorithms.html

    - Return node properties: dictionaries keyed by node label
    {REACTION: VALUE}

'''

import os
import sys
import networkx as nx


def progress(count, total):
    '''
    Show progress bar.
    Modified from https://gist.github.com/vladignatyev/06860ec2040cb497f0f3
    '''
    
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s\r' % (bar, percents, '%'))
    sys.stdout.flush()  
    
    
    
def make_folder(folder):
    if not os.path.exists(folder):
            os.makedirs(folder)
    
    

def create_directed_RG(out):
    '''
    Create a DIRECTED Reaction Graph
    '''
    edges_model = open(out + '/edge.list', 'rb')
    nodes_model = open(out + '/node.list', 'rb')
    nodes_model2 = nodes_model.read().splitlines()

    DG = nx.read_edgelist(edges_model, create_using=nx.DiGraph(), delimiter = '\t')
    DG.add_nodes_from(nodes_model2)
    DG.name = 'ReactionGraph'
    DG = nx.freeze(DG)
    ## Some stats about the graph
    f = open(out + '/stats.txt','w')
    f.write(nx.info(DG)+'\n')
    f.close()    
    return(DG)


def get_connected_components(out, DG):
    '''
    To obtain the connected components we need to create an INDIRECTED graph from 
    a DIRECTED graph. To transform it again to a DIRECTED graph we need to remove
    the duplicated links.     
    Returns a sorted list of DIRECTED connected components (Giant component first). 
    Writes: 
        - edgelist for each component --> read with: nx.read_edgelist(filename,create_using=nx.DiGraph())
        - nodelist for each component: isolated nodes have empty edgelist.
    
    '''   
    UDG = DG.to_undirected()
    print('\nCalculating number of connected components...')
    CC = list(sorted(nx.connected_component_subgraphs(UDG),key=len, reverse=True))

    if CC :
        count=0
        for mini in CC:
          miniD = mini.to_directed()  
          name = CC.index(mini)
          miniD.name = ('%03d' % name ) 
          
          # remove links that were not in the original files
          edges2remove = list(set(miniD.edges()) - set(DG.edges()))
          miniD.remove_edges_from(edges2remove)
          
          ## write edge list, node list in new folder
          newout = out+'/'+str(miniD.name)
          make_folder(newout)
          nx.write_edgelist(miniD, newout+'/edge.list', delimiter="\t", data=False)
          nodes_file = open(newout+'/node.list', 'w')
          for nodes in miniD.nodes():
              nodes_file.write(nodes+'\n')
          nodes_file.close()    
          # print progress
          count+=1
          progress(count, len(CC))

    print('\nNumber of connected components: '+str(len(CC)))

    return(CC)


def DG_indegree(out,DGc):
    '''
    Normalized indegree (values are normalized by /n-1) 
    '''
    make_folder(out+'/topology')
    measure = nx.in_degree_centrality(DGc)
    f = open(out+'/topology/indegree.list','w')
    f.write('REACTION\tINDEGREE\n')
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('indegree calculated')
    return(measure)

    
def DG_outdegree(out,DGc):
    '''
    Normalized outdegree (values are normalized by /n-1) 
    '''
    make_folder(out+'/topology')
    measure = nx.out_degree_centrality(DGc)
    f = open(out+'/topology/outdegree.list','w')
    f.write('REACTION\tOUTDEGREE\n')    
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('outdegree calculated')
    return(measure)
    
    
def DG_degree(out,DGc):
    '''
    Degree (values are normalized by /n-1) 
    '''
    make_folder(out+'/topology')
    measure = nx.degree_centrality(DGc)
    f = open(out+'/topology/degree.list','w')
    f.write('REACTION\tDEGREE\n')    
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('degree calculated')
    return(measure)

    
    
def DG_closeness(out,DGc):
    '''
    Closeness reciprocal of the sum of the shortest path distances from a node to all other nodes.
    If the graph is not completely connected, this algorithm computes the closeness centrality 
    for each connected part separately.
    '''
    make_folder(out+'/topology')
    measure = nx.closeness_centrality(DGc)                     
    f = open(out+'/topology/closeness.list','w')
    f.write('REACTION\tCLOSENESS\n')    
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('closeness calculated')
    return(measure)
    
    
def DG_betweenness(out,DGc):
    '''
    Betweeness
    '''
    make_folder(out+'/topology')
    measure = nx.betweenness_centrality(DGc)                     
    f = open(out+'/topology/betweenness.list','w')
    f.write('REACTION\tBETWEENESS\n')
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('betweenness calculated')
    return(measure)
    

def DG_eigen_left_ORI(out,DGc):
    '''
    Left Eigenvector centrality: For directed graphs this is "left" eigenvector centrality which 
    corresponds to the in-edges in the graph.     
    Using Numpy calculation to avoid non-convergence
    '''
    make_folder(out+'/topology')
    measure = nx.eigenvector_centrality_numpy(DGc)
    f = open(out+'/topology/eigen_left.list','w')
    f.write('REACTION\tEIGENLEFT\n')    
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('left eigenvector calculated')
    return(measure)



def DG_eigen_right_ORI(out,DGc):
    '''
    Right Eigenvector centrality: For out-edges eigenvector centrality first reverse the graph with G.reverse().
    Using Numpy calculation to avoid non-convergence
    '''
    make_folder(out+'/topology')
    DGcREV = DGc.reverse()
    measure = nx.eigenvector_centrality_numpy(DGcREV)
    f = open(out+'/topology/eigen_right.list','w')
    f.write('REACTION\tEIGENRIGHT\n')    
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('right eigenvector calculated')
    return(measure)




def DG_eigen_left(out,DGc):
    '''
    Left Eigenvector centrality: For directed graphs this is "left" eigenvector centrality which 
    corresponds to the in-edges in the graph.     
    Using Numpy calculation to avoid non-convergence
    '''
    warnings.filterwarnings("error")    
    make_folder(out+'/topology')
    try:
        measure = nx.eigenvector_centrality_numpy(DGc)
        f = open(out+'/topology/eigen_left.list','w')
        f.write('REACTION\tEIGENLEFT\n')    
        for k,v in measure.items():
            f.write(k+'\t'+str(v)+ '\n')
        f.close
        print('left eigenvector calculated')
    except RuntimeWarning:
        print('WARNING: L.Eigenvector calculation did not converge (score = NA)')
        measure = dict(zip(DGc.nodes(),['NA']*DGc.number_of_nodes() ) )     
    return(measure)

    
    
def DG_eigen_right(out,DGc):
    '''
    Right Eigenvector centrality: For out-edges eigenvector centrality first reverse the graph with G.reverse().
    Using Numpy calculation to avoid non-convergence
    '''
    warnings.filterwarnings("error")    
    make_folder(out+'/topology')
    DGcREV = DGc.reverse()
    try:
        measure = nx.eigenvector_centrality_numpy(DGcREV)
        f = open(out+'/topology/eigen_right.list','w')
        f.write('REACTION\tEIGENRIGHT\n')    
        for k,v in measure.items():
            f.write(k+'\t'+str(v)+ '\n')
        f.close
        print('right eigenvector calculated')
    except RuntimeWarning:
        print('WARNING: R.Eigenvector calculation did not converge (score = NA)')
        measure = dict(zip(DGcREV.nodes(),['NA']*DGcREV.number_of_nodes() ) )     

    return(measure)



def DG_ratio_io(out,DGc):
    '''
    Calculate the in/out-degree ratio by node. 
    Careful with ZeroDivisionError
    '''
    make_folder(out+'/topology')
    measure = {}
    for n in DGc.nodes():
          outd = DGc.in_degree(n)     # indegree not normalized
          ind = DGc.out_degree(n)     # outdegree not normalized
          if outd == 0:
            ratio = 0.0
          else:
            ratio = float(ind)/outd   # ratio in/outdegree
          measure[n] = ratio
    f = open(out+'/topology/ratio_io.list','w')
    f.write('REACTION\tRATIO_IO\n')
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('ratio in/outdegree calculated')
    return(measure)



def DG_source_sink(out,DGc):
    '''
    Compute Source/Sink (Input/Yeld) categories based on node degree:
      - Source (input): indegree = 0
      - Sink (yeld): outdegree = 0 
      - Intermediate (noinput_noyeld): outdegree > 0 & indegree > 0
      - Isolated (input_yeld): outdegree = 0 & indegree = 0 
      
    In a connected component there wouldn't be isolated nodes.    

    Modified from Ludovica Montanucci yeld_input function.
    '''
    make_folder(out+'/topology')
    measure = {}
    indegree = nx.in_degree_centrality(DGc)
    outdegree = nx.out_degree_centrality(DGc)
    for node in indegree.keys():
        if indegree[node] == 0.0 and outdegree[node] != 0.0:
            measure[node] = 'source' #'input'
        elif indegree[node] != 0.0 and outdegree[node] == 0.0:
            measure[node] = 'sink' #'yeld' 
        elif indegree[node] != 0.0 and outdegree[node] != 0.0:
            measure[node] = 'intermediate' #'noinput_noyeld' 
        elif indegree[node] == 0.0 and outdegree[node] == 0.0:
            measure[node] = 'isolated' #'input_yeld' 

    f = open(out+'/topology/source_sink.list','w')
    f.write('REACTION\tSOURCE_SINK\n')    
    for k,v in measure.items():
        f.write(k+'\t'+str(v)+ '\n')
    f.close
    print('source/sink calculated')
    return(measure)



def DG_get_num_pred_succ(dict_pred_succ):
  '''
  Calculate the number of predecessors or successors from a dictionary. 
  Cannot do it directly because of recursive loop in previous function.
  '''
  dictionary = dict_pred_succ
  for k in dictionary:
    dictionary[k] = len(dictionary[k])
  return dictionary



def DG_get_predecessors(DirG,node,all_predec):
  '''
  Networkx function only gets the immediate predecessor node of n. 
  Modified from Ludovica Montanucci.
  '''
  for n in DirG.predecessors_iter(node):
    all_predec=all_predec.union(DG_get_predecessors(DirG,n,all_predec))
    all_predec=all_predec.union(set([n]))
  return all_predec



def DG_get_successors(DirG,node,all_succ):
  '''
  Networkx function only gets the immediate successor node of n. 
  Modified from Ludovica Montanucci.
  '''
  for n in DirG.successors_iter(node):
    all_succ=all_succ.union(DG_get_successors(DirG,n,all_succ))
    all_succ=all_succ.union(set([n]))
  return all_succ



def DG_successors_predecessors(out,DGc):
    '''
    In a linear pathway, count:
        -Successors: number of nodes(reactions) that go after that node. 
        -Predecessors: number of nodes(reactions) that go before that node.
    
    A linear pathway is a pathway without feedback loops or reversible reactions.
    Remove self-feedback loops to avoid an infinite recursive loop.      
    Modified from Ludovica Montanucci.  
    '''
    make_folder(out+'/topology')
    measure_pred, measure_succ,num_measure_pred,num_measure_succ = {}, {}, {}, {}
    DirG = nx.DiGraph(DGc) # need to copy it to remove self-feedback loops
    DirG.remove_edges_from(DirG.selfloop_edges())  

    n_simple_cycles = nx.simple_cycles(DirG) 
    if next(n_simple_cycles,None):    # NOT directly linear (reversible rr, loop)
        print('WARNING: Graph has cycles (successors/predecessors = CYC)')
        #num_measure_pred = dict(zip(DirG.nodes(),['CYC']*DirG.number_of_nodes() ) ) 
        #num_measure_succ = dict(zip(DirG.nodes(),['CYC']*DirG.number_of_nodes() ) ) 
    else:
        for node in DirG.nodes():       # LINEAR GRAPH 
            measure_succ[node] = DG_get_successors(DirG,node, set([]) )
            measure_pred[node] = DG_get_predecessors(DirG,node, set([]) )
        num_measure_succ = DG_get_num_pred_succ(measure_succ)
        num_measure_pred = DG_get_num_pred_succ(measure_pred)
        print('successors/predecessors calculated')
    
    if len(num_measure_succ) > 0:    
        f = open(out+'/topology/sucessors.list','w')
        f.write('REACTION\tSUCCESORS\n')
        for k,v in num_measure_succ.items():
            f.write(k+'\t'+str(v)+ '\n')
        f.close
        f = open(out+'/topology/predecessors.list','w')
        f.write('REACTION\tPREDECESSORS\n')
        for k,v in num_measure_pred.items():
            f.write(k+'\t'+str(v)+ '\n')
        f.close
        return(num_measure_succ, num_measure_pred)


if __name__ == '__main__':

    ## Get arguments
    ifiles = sys.argv[1] # files inside reaction_graph
    output = sys.argv[2] # cComponents
    
    ## Create the main directed Reaction Graph
    DirRG = create_directed_RG(ifiles)

    ## Generate all connected components of the graph (undirected)
    graphs = get_connected_components(output, DirRG) 

    ## CALCULATE TOPOLOGICAL MEASURES - DIRECTED graph. ADD functions HERE
    if graphs:
        ccomp_list = sorted(os.listdir(output))
        for item in ccomp_list:
            print('\nComponent: '+item)
            newout = output+'/'+str(item)
            comp = create_directed_RG(newout)
            if comp.number_of_edges() > 0:
                DG_indegree(newout,comp)
                DG_outdegree(newout,comp)
                DG_ratio_io(newout,comp)
                DG_closeness(newout,comp)
                DG_betweenness(newout,comp)
                DG_source_sink(newout,comp) # source-sink-intermediate
                DG_successors_predecessors(newout,comp) # returns 2 dictionaries

                ## some components do not converge
                #if comp.number_of_nodes() > 2:
                #    DG_eigen_right(newout,comp)
                #    DG_eigen_left(newout,comp)
                #else:
                #    print('Eigenvector calculation needs nodes > 2')

            else:
                print('No edges found!')
    else:
        print('No connected components found')
        

