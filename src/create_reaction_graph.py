#!/usr/bin/env python

'''

Create a reaction graph from a metabolic model in MATLAB format.
The rules are: 

    - If the products of a node are the reactants of any other node, 
    create a directed link between them.
    - Removing currency metabolites.
    - Remove generic biomass reaction
    - Allowing self-loops. 
    - If reaction is reversible, check reactants signs (+/-) to identify phyisiological direction. 
    - Returns a series of files with information about the graph:
        - List of nodes (REACTIONS)
        - List of edges (DIRECTED)
        - List of genes (EntrezGene IDs)
        - List of subsystems (METABOLIC PATHWAYS)
        - Link between genes IDs and reactions
        
'''
import os
import sys
import getopt
import re
import cobra
import pandas as pd


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
    

def make_pwy_file(out,mod):
    '''
    Write a file with all metabolic subsystems ( = pathways) in the model. 
    '''  
    f = open(out + '/subsystems.list', 'w')
    list_of_pathways = []

    for r in mod.reactions: 
      pwy_r = r.subsystem
      if pwy_r.startswith('[array'): ## Recon3D has a different format that Recon2
            pwy_r = r.subsystem.split("\'")[1]
            list_of_pathways.append(pwy_r) 
      else:
            list_of_pathways.append(pwy_r) 

    for item in set(list_of_pathways):
      f.write("%s\n" % item)
    f.close()
    return set(list_of_pathways)
    

    
def make_geneReaction_file(out,mod,geneID='EntrezGene'):
    '''
    Write a file linking reactions and genes.
    If the gene participates in a reaction, it doesn't matter if it's the only one, 
    part of a subunit or an isoenzyme. 
   
    GeneID:
        - Entrez --> genes are coded with entrezGene transcripts ids (27349.1), keep only gene unique ids (27349)
        - Symbol --> unique ids
    Returns a pandas data.frame: GENE REACTION

    '''
    RXNsByGENE = {}
    df_g, df_r = [], []
    for g in mod.genes:
      g_id = g.id   # 27349.1 
      all_rr = []
      for rr in list(g.reactions):
        all_rr.append(rr.id)
      RXNsByGENE[g_id] = all_rr
    
    if geneID == 'EntrezGene':
        all_genes = set([ x.split('.')[0] for x in RXNsByGENE.keys() ] )
        for uniq_g in list(all_genes):
          values_rr = [ v for k,v in RXNsByGENE.items() if k.startswith(str(uniq_g+'.'))] 
          uniq_rr = set([item for sublist  in values_rr for item in sublist] )  
          for item in uniq_rr:
            df_r.append(item)
            df_g.append(uniq_g)
        genereactions_file = pd.DataFrame({'GENE': df_g, 'REACTION': df_r})
        genereactions_file.to_csv(out +'/geneReactions.list', sep='\t',index=False)
    elif geneID == 'Symbol':
        all_genes = set(RXNsByGENE.keys())
        for uniq_g in list(all_genes):
          values_rr = [ v for k,v in RXNsByGENE.items() if k == uniq_g] 
          uniq_rr = set([item for sublist in values_rr for item in sublist] )  
          for item in uniq_rr:
            df_r.append(item)
            df_g.append(uniq_g)
        genereactions_file = pd.DataFrame({'GENE': df_g, 'REACTION': df_r})
        genereactions_file.to_csv(out +'/geneReactions.list', sep='\t',index=False)
    else: 
        print('geneID = EntrezGene | Symbol | define your own!')

    return(genereactions_file)

    

def remove_currency_meta(mod, reaction, metas,remove_currency):
  '''
  Do not take into account currency metabolites with the highest degree: 
  16 metabolites in 8 compartments (128 metabolites)
  NOTE: compartments might change! Extract them directly!
  
  '''
  
  currency =  ["adp", "atp", "co2", "o2", "h2o", "h2o2", "h", "k", "na1", "nad", "nadh", "nadp", "nadph", "nh4", "pi", "ppi"]
  compartments = ["[c]","[e]", "[l]", "[m]", "[x]", "[r]", "[g]", "[n]"]
  currency_comp = []
  if remove_currency == True: 
    for c in currency:
      cc = ';'.join('%s%s' % (c,cm)  for cm in compartments)
      cc = cc.split(';')
      currency_comp.extend(cc)
    for m in metas:
      if str(m) in currency_comp:
        reaction.pop(mod.metabolites.get_by_id(str(m)))
  else:
    pass
  return(reaction)



def extract_metabolites_ordered(mod, reaction, remove_currency=True):
  '''
  For every reaction extract the metabolites (reactants, products)
  In reversible reactions --> metabolites -/+ sign indicates the physiological direction
  
  '''
  products_interch, reactants_interch = [],[]
  rr = mod.reactions.get_by_id(reaction)
  #remove currency metabolites from the reaction
  cleaned_reaction = remove_currency_meta(mod, rr,rr.metabolites,remove_currency)
  if cleaned_reaction.reversibility == True:
    for sign in cleaned_reaction.metabolites:
      if cleaned_reaction.metabolites[sign] == 1.0: # product
        products_interch.append(sign)
      elif cleaned_reaction.metabolites[sign] == -1.0: # reactant
        reactants_interch.append(sign)
    react_cleaned = reactants_interch
    prod_cleaned = products_interch
  else:
    react_cleaned = cleaned_reaction.reactants
    prod_cleaned = cleaned_reaction.products
  return(react_cleaned,prod_cleaned)
    

def get_all_nodes_ids(all_model_reactions):
  '''
  Extract all nodes by reactions ids
  '''
  list_of_nodesIDS = []
  for node in all_model_reactions:
    list_of_nodesIDS.append(node.id)
  return(list_of_nodesIDS)



def make_links(mod, all_model_reactions, remove_currency=True):
  '''
  If the products of node are the reactants of other nodes --> create a link between them
  Allow for self-loops 

  '''
  list_of_self_edges = []
  list_of_edges = []
  count=0
  print('\nCalculating links...')  
  for node1 in all_model_reactions:
    #Extract reactants and products for every reaction (node)  
    react_cleaned,prod_cleaned = extract_metabolites_ordered(mod, node1.id, remove_currency)
    # Make a self-loop
    if any(map(lambda y: y in prod_cleaned, react_cleaned)):
      list_of_self_edges.append(node1.id+'\t'+node1.id)
    # Get IDs of all other nodes
    other_nodes = get_all_nodes_ids(mod.reactions)      
    other_nodes.remove(node1.id)
    for nodes_by_id in other_nodes:
      node2 = mod.reactions.get_by_id(nodes_by_id)
      #Extract reactants and products for every reaction (node)  
      react_cleaned2,prod_cleaned2 = extract_metabolites_ordered(mod, node2.id, remove_currency)
      # If the products of node1 are the reactants of other nodes --> create a link between them
      if any(map(lambda x: x in prod_cleaned, react_cleaned2)):
        list_of_edges.append(node1.id+'\t'+node2.id)
    # print progress
    count+=1
    progress(count, len(all_model_reactions))

  all_edges = list_of_edges + list_of_self_edges
  all_edges = set(all_edges)  
  return(all_edges)



def make_edge_file(out,mod):
    '''
    Write a file with all edges between nodes, edges are directed. 
    '''
    links_model = make_links(mod, mod.reactions, remove_currency = True)

    f = open(out +'/edge.list', 'w')
    for item in links_model:
        f.write("%s\n" % item)
    f.close()
    return(links_model)



def make_node_file(out,mod):
    '''
    Write a file with all nodes, where nodes are the reactions
    '''
    f = open(out +'/node.list', 'w')
    for nodes in mod.reactions:
      f.write("%s\n" % nodes.id)
    f.close()
    return(mod.reactions)



def make_gene_file(out,mod, geneID='EntreGene'):
    '''
    Write a file with genes in the model.
    GeneID:
        - Entrez --> genes are coded with entrezGene transcripts ids (27349.1), keep only gene unique ids (27349)
        - Symbol --> unique ids
    '''
    
    
    if geneID == 'EntrezGene':
        all_genes_dup = []
        for g in mod.genes:
            all_genes_dup.append(g.id)      
        all_genes = set([ x.split('.')[0] for x in all_genes_dup] )
    elif geneID == 'Symbol':
        all_genes = []
        for g in mod.genes:
            all_genes.append(g.id)      
    else: 
        print('geneID = EntrezGene | Symbol | define your own!')
  
    f = open(out +'/gene.list', 'w')
    for g in all_genes:
      f.write("%s\n" % g)
    f.close()
    return(all_genes)



if __name__ == '__main__':

    ## Get arguments 
    output = sys.argv[1]
    matfile = sys.argv[2]
       
    model = cobra.io.load_matlab_model(matfile)

    ## NOTE: removing generic biomass reaction from model
    biomass= [r.id for r in model.reactions if re.search('biomass', r.id,re.IGNORECASE)]
    print('\nRemoving reactions: '+str(biomass))
    if biomass:
        model.remove_reactions(biomass)

    ## Create file with nodes --> REACTION
    nodesModel = make_node_file(output, model)
    print('\nNumber of nodes: '+str(len(nodesModel)))

    ## Create file with edges (directed)
    edgesModel = make_edge_file(output, model)
    print('\nNumber of links: '+str(len(edgesModel)))

    ## Create file with subsystems (Pathways)
    subsystemsModel = make_pwy_file(output, model)
    print('\nNumber of subsystems: '+ str(len(subsystemsModel)))

    ## Create file with genes (geneID: EntrezGene | Symbol )
    genesModel = make_gene_file(output, model, geneID = 'EntrezGene')
    print('\nNumber of genes:' +str(len(genesModel)))

    ## Create file linking reactions and genes (geneID: EntrezGene | Symbol )
    geneReaction = make_geneReaction_file(output, model, geneID = 'EntrezGene')

