import pulp
import json
#import pandas as pd
from tqdm import tqdm
from rdkit import Chem
import numpy as np
import time
import streamlit as st

@st.cache_data
def load_metab_df():
    metab_df = pd.read_csv("./metanetx/data_final/metanetx_metab_db_noduplicates.csv" , index_col = "Unnamed: 0")
    return metab_df
    
@st.cache_data
def load_files():
    
    all_sij = json.load(open("./all_sij_90.json"))
    all_rij_NEW = json.load(open("./all_rij_90.json"))
    rxn_classify = json.load(open("./../../rxn_classify.json"))
    #bio_rxns_front = json.load(open("./../updated_data/new_data/bio_rxns_front.json"))
    rev_pair = json.load(open("./rev_pair_90.json"))
    new_all_mol_ids = list(all_rij_NEW.keys())
    new_all_rxn_ids = list(all_sij.keys())
    bio_chem_smiles_ids_dict = json.load(open('./../../bio_chem_smiles_ids_dict_updated.json'))
    
    return new_all_mol_ids, new_all_rxn_ids, rxn_classify, all_sij,all_rij_NEW,bio_chem_smiles_ids_dict,rev_pair


def add_exchange(exchange_mets,bio_chem_smiles_ids_dict,new_all_rxn_ids,rxn_classify,all_rij_NEW):
    exchange_rxns = ["EX_"+i for i in exchange_mets]
    all_rij_NEW[exchange_mets[0]][exchange_rxns[0]] = 1
    all_rij_NEW[exchange_mets[1]][exchange_rxns[0]] = -1
    for i,ids in enumerate(exchange_rxns):
        new_all_rxn_ids.append(ids)
        rxn_classify[ids] = 2
    
    return new_all_rxn_ids,rxn_classify,all_rij_NEW,exchange_rxns
def formulation(new_all_mol_ids,new_all_rxn_ids,all_rij_NEW,rxn_classify,exchange_rxns,soln_file,tot_transitions,rev_pair):
    M_low = 1
    M_high = 1
    eps = 0
    start = time.time()
    #Variables
    #y_rxn = pulp.LpVariable.dicts("y_rxn", new_all_rxn_ids, cat="Integer")        
    #y_rxn = pulp.LpVariable.dicts("y_rxn", new_all_rxn_ids, lowBound=M_low, upBound=M_high, cat="Binary")
    y_rxn = pulp.LpVariable.dicts("y_rxn", new_all_rxn_ids, lowBound=0, upBound=1, cat="Binary")
    y_B = pulp.LpVariable.dicts("y_B", new_all_mol_ids, lowBound=0, upBound=1, cat="Binary")
    y_B2 = pulp.LpVariable.dicts("y_B2", new_all_mol_ids, lowBound=0, upBound=1, cat="Binary")
    y_C = pulp.LpVariable.dicts("y_C", new_all_mol_ids, lowBound=0, upBound=1, cat="Binary")
    y_C2 = pulp.LpVariable.dicts("y_C2", new_all_mol_ids, lowBound=0, upBound=1, cat="Binary")
    y_BC = pulp.LpVariable.dicts("y_BC", new_all_mol_ids, lowBound=0, upBound=1, cat="Binary")
    y_CB = pulp.LpVariable.dicts("y_CB", new_all_mol_ids, lowBound=0, upBound=1, cat="Binary")
    #Objective function
    lp_prob = pulp.LpProblem("chemoenzymatic", pulp.LpMinimize)
    lp_prob += pulp.lpSum([y_rxn[j] for j in new_all_rxn_ids if rxn_classify[j]==1])+pulp.lpSum([y_BC[i]+y_CB[i] for i in new_all_mol_ids]) + pulp.lpSum([0.5*y_rxn[j] for j in new_all_rxn_ids if rxn_classify[j]==0])
    #Constraints on y_rxn
    '''
    for rids in new_all_rxn_ids:
        if rxn_classify[rids]==0:
            lp_prob += y_rxn[rids]>=-1, "y_rxn_"+rids+"_lower_bound"
            lp_prob += y_rxn[rids]<=1, "y_rxn_"+rids+"_upper_bound"
        else:
            lp_prob += y_rxn[rids]>=0, "y_rxn_"+rids+"_lower_bound"
            lp_prob += y_rxn[rids]<=1, "y_rxn_"+rids+"_upper_bound"
    '''
    #Constraints
    
    for rids in rev_pair:
        all_combs = rev_pair[rids]
        for combs in all_combs:
            lp_prob+= y_rxn[rids]+y_rxn[combs]<=1, "reversible_"+rids+"_"+combs
        
    #'MNXR150321_1_rev', 'MNXR188567_rev'
    #'MNXM1103718' -> CTP, 'MNXM1102191' -> CDP, 'MNXM411' -> dCDP,'MNXM360->dCTP
    #'MNXM1102128' -> UDP, 'MNXM1101474' -> UTP, 'MNXM728294' -> AMP, 'MNXM40333'- ADP, 'MNXM1103302' -> CMP
    #{'MNXM1105937': FAD , 'MNXM1105762':FADH2}
    # ->  
    #cofac_list = ['MNXM1103718','MNXM1102191','MNXM411','MNXM1102128','MNXM1101474','MNXM728294','MNXM40333', 'MNXM1103302','MNXM1105937','MNXM1105762']
    cofac_met = ['MNXM40333','MNXM1103718','MNXM1101474','MNXM1103428','MNXM1102128','MNXM1103285','MNXM1102191','MNXM1104559','MNXM728294','MNXM1103302','MNXM1104823','MNXM1101285','MNXM1101868','MNXM728062','MNXM1102072','MNXM1105762','MNXM411','MNXM1105937']

    cofac_more = ['MNXM152','MNXM436', 'MNXM1102167', 'MNXM735047', 'MNXM1094084', 'MNXM1103553', 'MNXM191', 'MNXM232', 'MNXM1108515', 'MNXM13204', 'MNXM731831', 'MNXM723089', 'MNXM113802', 'MNXM123', 'MNXM736654', 'MNXM8529','MNXM740811']
    for ids in cofac_met:
        lp_prob+=y_C[ids]==0, "cofac_"+ids
    for ids in cofac_more:
        lp_prob+=y_C[ids]==0, "cofac_"+ids
    
    #lp_prob+=y_rxn['MNXR165099_3_rev']==0, "oxoglutarate_useless"
    #lp_prob+=y_rxn['MNXR188567_rev']==0, "infeasible_DG"
    
    #lp_prob+=y_C['MNXM1103302']==0, "CMP"
    
    
    '''
    for rids in new_all_rxn_ids:
        if "M" in rids:
            rids_rev = rids + "_rev"
            if rids_rev in new_all_rxn_ids:
                lp_prob += y_rxn[rids]+y_rxn[rids_rev]<=1, "reversible_"+rids
    '''
    
    for ids in new_all_mol_ids:
        lp_prob += pulp.lpSum([all_rij_NEW[ids][j]*y_rxn[j] for j in list(all_rij_NEW[ids].keys())]) == 0, "molecule_balance_"+ids
    '''
    for j in new_all_rxn_ids:
        lp_prob += y_rxn[j] >= y_rxn[j]* M_low, "cons1_" + j
        lp_prob += y_rxn[j] <= y_rxn[j]* M_high, "cons2_" + j
   '''
    #0 for B, 1 for C
    for ids in new_all_mol_ids:

        #C to B transition

        lp_prob += pulp.lpSum([-all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]<0 and rxn_classify[j]==1]) <= y_B[ids]*M_high, "cons3_" + ids
        lp_prob += pulp.lpSum([-all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]<0 and rxn_classify[j]==1]) >= y_B[ids]*M_low, "cons4_" + ids

        lp_prob += pulp.lpSum([all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]>0 and rxn_classify[j]==0]) <= y_C[ids]*M_high, "cons5_" + ids
        lp_prob += pulp.lpSum([all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]>0 and rxn_classify[j]==0]) >= y_C[ids]*M_low, "cons6_" + ids


        lp_prob += y_BC[ids] <= y_B[ids], "cons11_" + ids
        lp_prob += y_BC[ids] <= y_C[ids], "cons12_" + ids
        lp_prob += y_BC[ids] >= y_B[ids] + y_C[ids] - 1, "cons13_" + ids

        #B to C transition

        lp_prob += pulp.lpSum([all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]>0 and rxn_classify[j]==1]) <= y_B2[ids]*M_high, "cons7_" + ids
        lp_prob += pulp.lpSum([all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]>0 and rxn_classify[j]==1]) >= y_B2[ids]*M_low, "cons8_" + ids

        lp_prob += pulp.lpSum([-all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]<0 and rxn_classify[j]==0]) <= y_C2[ids]*M_high, "cons9_" + ids
        lp_prob += pulp.lpSum([-all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]<0 and rxn_classify[j]==0]) >= y_C2[ids]*M_low, "cons10_" + ids


        lp_prob += y_CB[ids] <= y_B2[ids], "cons14_" + ids
        lp_prob += y_CB[ids] <= y_C2[ids], "cons15_" + ids
        lp_prob += y_CB[ids] >= y_B2[ids] + y_C2[ids] - 1, "cons16_" + ids
        
        # Every mol should be produced once & consumed once
        
        #lp_prob += pulp.lpSum([all_rij_NEW[ids][j]*y_rxn[j] for j in all_rij_NEW[ids] if all_rij_NEW[ids][j]!=0) <= 2, "cons17_" + ids
    
        
    #lp_prob += y_rxn[exchange_rxns[0]]==1, "input_flux"
    #lp_prob += y_rxn[exchange_rxns[1]]==1, "output_flux"
    lp_prob += y_rxn[exchange_rxns[0]]==1, "exchange_input"
    lp_prob += y_rxn[exchange_rxns[1]]==1, "exchange_output"
    #lp_prob += y_rxn['MNXR153413']==1, "making mep"
    #lp_prob+= pulp.lpSum([y_rxn[j] for j in new_all_rxn_ids if rxn_classify[j]!=2]) >= 12 ,"minimum_steps"
    lp_prob+= pulp.lpSum([y_rxn[j] for j in new_all_rxn_ids if rxn_classify[j]!=2]) <= 30, "maximum_steps"
    #lp_prob += pulp.lpSum([y_BC[ids]+y_CB[ids] for ids in new_all_mol_ids])==0, "transitions"
    itr = 0
    for num_transitions in range(tot_transitions):
        integer_cut_rules = []
        integer_cut_vals = []
        #lp_prob.constraints['transitions'].changeRHS(num_transitions)
        itr+=1
        #for itr in range(5):
        #path = "./log_files_FPP_2/"+str(num_transitions)+"_"+str(itr)+".log"
        #pulp_solver = pulp.CPLEX_CMD(path=None,keepFiles=0, mip=1, msg=1, logPath=path)
        #pulp_solver = pulp.CPLEX_CMD(path=None,keepFiles=0, mip=1, msg=1, logPath=path)
        pulp_solver = pulp.CPLEX_CMD(path=None,keepFiles=0, mip=1, msg=1, timeLimit=3000)
        pulp.LpSolverDefault.msg = 1
        print("---------going to the solvers--------")
        lp_prob.solve(pulp_solver)
        print("---------DId not get any error--------")
        end = time.time() 
        print("-----------------------------")
        with open(soln_file,'a') as f:
            f.write('\n------------------------------------')
            f.write('\nTime = {}'.format(end-start))
            f.write("\nIteration num = {}".format(itr))
            f.write("\nStatus = {}\n".format(pulp.LpStatus[lp_prob.status]))
            f.write("\nTransitions = {}".format(num_transitions))
            f.write("\n")
        #print("Status:", pulp.LpStatus[lp_prob.status])
        #lp_prob.writeMPS("./log_files_tolyl_alt/mps/"+str(num_transitions)+"_"+str(itr)+".mps")
        #lp_prob.writeLP("./log_files_tolyl_alt/lp/"+str(num_transitions)+"_"+str(itr)+".lp")
        if pulp.LpStatus[lp_prob.status] != 'Optimal':
            break
        else:
            #pulp.value(my_lp_problem.objective)
            obj = pulp.value(lp_prob.objective)
            with open(soln_file,'a') as f:
                f.write("Objective value = {}\n".format(obj))
                f.write('y_rxn_values')
                f.write('\n')
            for rid in new_all_rxn_ids:
                if y_rxn[rid].varValue != 0:
                    integer_cut_rules.append(rid)
                    integer_cut_vals.append(y_rxn[rid].varValue)
                    #print(rid + ',' + str(y_rxn[rid].varValue)+" y_rxn =")
                    #print(str(y_rxn[rid].varValue))
                    #print('\n')
                    with open(soln_file,'a') as f:
                        f.write(rid + ',' + str(y_rxn[rid].varValue)+" y_rxn =")
                        f.write(str(y_rxn[rid].varValue))
                        f.write('\n')
                    #if obj

            with open(soln_file,'a') as f:
                f.write('\n')
                f.write('y_CB_values')
                f.write('\n')
            for ids in new_all_mol_ids:
                if y_CB[ids].varValue != 0:
                    with open(soln_file,'a') as f:
                        f.write(str(y_CB[ids].name) + ',' + str(y_CB[ids].varValue))
                        f.write('\n')   
            with open(soln_file,'a') as f:
                f.write('\n')
                f.write('y_BC_values')
                f.write('\n')
            for ids in new_all_mol_ids:
                if y_BC[ids].varValue != 0:
                    with open(soln_file,'a') as f:
                        f.write(str(y_BC[ids].name)+ ',' + str(y_BC[ids].varValue))
                        f.write('\n')
            '''
            with open(soln_file,'a') as f:
                f.write('\n')
                f.write('y_C_values')
                f.write('\n')           
            for ids in new_all_mol_ids:
                if y_C[ids].varValue !=0:
                    with open(soln_file,'a') as f:
                        f.write(str(y_C[ids].name)+ ',' + str(y_C[ids].varValue))
                        f.write('\n')
            with open(soln_file,'a') as f:
                f.write('\n')
                f.write('y_C2_values')
                f.write('\n')           
            for ids in new_all_mol_ids:
                if y_C2[ids].varValue !=0:
                    with open(soln_file,'a') as f:
                        f.write(str(y_C2[ids].name)+ ',' + str(y_C2[ids].varValue))
                        f.write('\n')
            '''       

            with open(soln_file,'a') as f:
                f.write('\n')
                f.write('Integer cut rules')
                f.write('{}\n\n'.format(integer_cut_rules))
                #f.write('{}\n\n'.format(integer_cut_vals))

            length = len(integer_cut_rules) - 1
            total_vals = sum(integer_cut_vals) - 1
            lp_prob += (
                pulp.lpSum([y_rxn[r] for r in integer_cut_rules]) <= total_vals,
                "integer_cut_" + str(itr)+"_trans_"+str(num_transitions),
            )
    
    return



def main():
    st.image('./../MinChemBio_header.png', use_column_width=True)
    st.subheader('Primary reactant & Primary product (Use MetaNetX IDs or My Chemsitry IDs)')
    #reactant = st.text_input('reactant', value='MNXM23')
    product = st.text_input(
            'product', value='CHEM01229150')

    
    
    if st.button("Search"):
        new_all_mol_ids, new_all_rxn_ids, rxn_classify, all_sij,all_rij_NEW,bio_chem_smiles_ids_dict,rev_pair = load_files()
        exchange_mets = ["MNXM23","CHEM01229150"]
        if product in all_rij_NEW:
            new_all_rxn_ids,rxn_classify,all_rij_NEW,exchange_rxns = add_exchange(exchange_mets,bio_chem_smiles_ids_dict,new_all_rxn_ids,rxn_classify,all_rij_NEW)
            st.write('Product is present in a reaction')
            metab_df = load_metab_df()
            metab_df_name = metab_df['Name'].to_dict()
            cid_fda_dict = json.load(open('./../cid_fda_dict.json'))
            # if session_state.button_search:
            st.subheader('Calculating chemoenzymatic pathway')
            st.write(reactant+" => "+ product)
            soln_file = "./results/"+cid_fda_dict[product]+"_"+metab_df_name[reactant]+"_.txt"
            tot_transitions = 100
            formulation(new_all_mol_ids,new_all_rxn_ids,all_rij_NEW,rxn_classify,exchange_rxns,soln_file,tot_transitions,rev_pair)
        else:
            st.write("PRODUCT NOT IN ANY REACTION")
            

if __name__ == '__main__':
    main()


