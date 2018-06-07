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
import itertools



def make_node_file(out,mod):
    '''
    Write a file with all nodes, where nodes are the reactions
    '''
    f = open(out +'/node.list', 'w')
    for nodes in mod.reactions:
      f.write("%s\n" % nodes.id)
    f.close()
    return(mod.reactions)
    
    
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
    
    

def make_gene_file(out,mod):
    '''
    Write a file with genes in the model.
    GeneID:
        - Entrez Gene --> genes are coded with entrezGene transcripts ids (27349.1), keep only gene unique ids (27349)
        - Ensembl Stable Gene ID --> unique ids
        - Other --> unique ids
    '''
    test_gene = mod.genes[0]
    
    if test_gene.id[0].isdigit():        
        all_genes_dup = []
        for g in mod.genes:
            all_genes_dup.append(g.id)      
        all_genes = set([ x.split('.')[0] for x in all_genes_dup] )
        print('\nGene ID seems EntrezGene')
    elif test_gene.id.startswith('ENSG'):
        all_genes = []
        for g in mod.genes:
            all_genes.append(g.id)  
        print('\nGene ID seems Ensembl Stable Gene ID')
    else:
        all_genes = []
        for g in mod.genes:
            all_genes.append(g.id)      
        print('\nGene ID is not detected as EntrezGene or Ensembl, creating file anyway')
  
    f = open(out +'/gene.list', 'w')
    for g in all_genes:
      f.write("%s\n" % g)
    f.close()
    return(all_genes)

    
def make_geneReaction_file(out,mod):
    '''
    Write a file linking reactions and genes.
    If the gene participates in a reaction, it doesn't matter if it's the only one, 
    part of a subunit or an isoenzyme. 
   
    GeneID:
        - Entrez Gene --> genes are coded with entrezGene transcripts ids (27349.1), keep only gene unique ids (27349)
        - Ensembl Stable Gene ID --> unique ids
        - Other --> unique ids
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
    
    test_gene = mod.genes[0]
    
    if test_gene.id[0].isdigit(): 
        all_genes = set([ x.split('.')[0] for x in RXNsByGENE.keys() ] )
        for uniq_g in list(all_genes):
          values_rr = [ v for k,v in RXNsByGENE.items() if k.startswith(str(uniq_g+'.'))] 
          uniq_rr = set([item for sublist  in values_rr for item in sublist] )  
          for item in uniq_rr:
            df_r.append(item)
            df_g.append(uniq_g)
        genereactions_file = pd.DataFrame({'GENE': df_g, 'REACTION': df_r})
        genereactions_file.to_csv(out +'/geneReactions.list', sep='\t',index=False)
    else:
        all_genes = set(RXNsByGENE.keys())
        for uniq_g in list(all_genes):
          values_rr = [ v for k,v in RXNsByGENE.items() if k == uniq_g] 
          uniq_rr = set([item for sublist in values_rr for item in sublist] )  
          for item in uniq_rr:
            df_r.append(item)
            df_g.append(uniq_g)
        genereactions_file = pd.DataFrame({'GENE': df_g, 'REACTION': df_r})
        genereactions_file.to_csv(out +'/geneReactions.list', sep='\t',index=False)

    return(genereactions_file)


def remove_currency_meta(mod, reaction, metas):
    '''
    Do not take into account currency metabolites with the highest degree: 
    16 metabolites in 8 compartments (128 metabolites)
    NOTE: compartments might change!
    '''
    #currency =  ["adp", "atp", "co2", "o2", "h2o", "h2o2", "h", "k", "na1", "nad", "nadh", "nadp", "nadph", "nh4", "pi", "ppi"]
    #compartments = ["[c]","[e]", "[l]", "[m]", "[x]", "[r]", "[g]", "[n]"]
    #currency_comp = []
    #for c in currency:
    #    cc = ';'.join('%s%s' % (c,cm)  for cm in compartments)
    #    cc = cc.split(';')
    #    currency_comp.extend(cc)
    currency_comp = ['adp[c]', 'adp[e]', 'adp[l]', 'adp[m]', 'adp[x]', 'adp[r]', 'adp[g]', 'adp[n]', 'atp[c]', 'atp[e]', 'atp[l]', 'atp[m]', 'atp[x]', 'atp[r]', 'atp[g]', 'atp[n]', 'co2[c]', 'co2[e]', 'co2[l]', 'co2[m]', 'co2[x]', 'co2[r]', 'co2[g]', 'co2[n]', 'o2[c]', 'o2[e]', 'o2[l]', 'o2[m]', 'o2[x]', 'o2[r]', 'o2[g]', 'o2[n]', 'h2o[c]', 'h2o[e]', 'h2o[l]', 'h2o[m]', 'h2o[x]', 'h2o[r]', 'h2o[g]', 'h2o[n]', 'h2o2[c]', 'h2o2[e]', 'h2o2[l]', 'h2o2[m]', 'h2o2[x]', 'h2o2[r]', 'h2o2[g]', 'h2o2[n]', 'h[c]', 'h[e]', 'h[l]', 'h[m]', 'h[x]', 'h[r]', 'h[g]', 'h[n]', 'k[c]', 'k[e]', 'k[l]', 'k[m]', 'k[x]', 'k[r]', 'k[g]', 'k[n]', 'na1[c]', 'na1[e]', 'na1[l]', 'na1[m]', 'na1[x]', 'na1[r]', 'na1[g]', 'na1[n]', 'nad[c]', 'nad[e]', 'nad[l]', 'nad[m]', 'nad[x]', 'nad[r]', 'nad[g]', 'nad[n]', 'nadh[c]', 'nadh[e]', 'nadh[l]', 'nadh[m]', 'nadh[x]', 'nadh[r]', 'nadh[g]', 'nadh[n]', 'nadp[c]', 'nadp[e]', 'nadp[l]', 'nadp[m]', 'nadp[x]', 'nadp[r]', 'nadp[g]', 'nadp[n]', 'nadph[c]', 'nadph[e]', 'nadph[l]', 'nadph[m]', 'nadph[x]', 'nadph[r]', 'nadph[g]', 'nadph[n]', 'nh4[c]', 'nh4[e]', 'nh4[l]', 'nh4[m]', 'nh4[x]', 'nh4[r]', 'nh4[g]', 'nh4[n]', 'pi[c]', 'pi[e]', 'pi[l]', 'pi[m]', 'pi[x]', 'pi[r]', 'pi[g]', 'pi[n]', 'ppi[c]', 'ppi[e]', 'ppi[l]', 'ppi[m]', 'ppi[x]', 'ppi[r]', 'ppi[g]', 'ppi[n]']

    [reaction.pop(mod.metabolites.get_by_id(str(m))) for m in metas if str(m) in currency_comp]

    return(reaction)


def extract_metabolites_ordered(mod, reaction, remove_currency):
    '''
    For every reaction extract the metabolites (reactants, products)
    In reversible reactions --> metabolites -/+ sign indicates the physiological direction

    '''
    products_interch, reactants_interch = [],[]
    rr = mod.reactions.get_by_id(reaction)
    #remove currency metabolites from the reaction?
    if remove_currency == True: 
        cleaned_reaction = remove_currency_meta(mod, rr,rr.metabolites)
    else:
        cleaned_reaction = rr
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



def make_links(node1, mod, remove_currency):
    '''
    If any of the products of a given node are the reactants of another node --> create a link between them
    Allows for self-loops: if the products of a node are also its reactants.
    Keep unique edges: if more than one product is a reactant count only 1 link.

    Extract reactants and products for a node: 
    react_cleaned,prod_cleaned = extract_metabolites_ordered(mod, node1.id, remove_currency)
    reactants = extract_metabolites_ordered(mod, node1.id, remove_currency)[0]
    products = extract_metabolites_ordered(mod, node1.id, remove_currency)[1]

    '''

    #Extract products for a node1   
    prod_cleaned = extract_metabolites_ordered(mod, node1.id, remove_currency)[1]
    
    #If any of the products of node1 are reactants of node make a link 
    list_of_edges = [
        node1.id+'\t'+node2.id
        for node2 in mod.reactions
        if any(map(lambda x: x in prod_cleaned, extract_metabolites_ordered(mod, node2.id, remove_currency)[0]))
    ]    

    all_edges_node = set(list_of_edges)  
    return(all_edges_node)



def make_edge_file(out,mod, rm_currency):
    '''
    Write a file with all edges between nodes, edges are directed. 
    With or without removing curreny metabolites to compare graphs.
    '''
    print('\nCalculating links...')  

    all_links_model_nodes = [ 
        make_links(nodes,mod,remove_currency=rm_currency) 
        for nodes in mod.reactions  
    ]
    
    all_links_model = list(itertools.chain(*all_links_model_nodes))

    if rm_currency is True:
        f = open(out +'/edge.list', 'w')
    else:
        f = open(out +'/edge_withCurrency.list', 'w')
    for item in all_links_model:
        f.write("%s\n" % item)
    f.close()
    return(all_links_model)




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

    ## Create file with edges (directed) --> keep currency metabolites
    edgesModelwCurrency = make_edge_file(output, model, rm_currency = False)
    print('\nNumber of links (with currency metabolites): '+str(len(edgesModelwCurrency)))
 
    ## Create file with edges (directed) --> remove currency metabolites
    edgesModel = make_edge_file(output, model, rm_currency = True)
    print('\nNumber of links (no currency metabolites): '+str(len(edgesModel)))
    
    ## Create file with subsystems (Pathways)
    subsystemsModel = make_pwy_file(output, model)
    print('\nNumber of subsystems: '+ str(len(subsystemsModel)))

    ## Create file with genes (geneID: EntrezGene | Ensembl )
    genesModel = make_gene_file(output, model)
    print('\nNumber of genes:' +str(len(genesModel)))

    ## Create file linking reactions and genes (geneID: EntrezGene | Ensembl )
    geneReaction = make_geneReaction_file(output, model)

