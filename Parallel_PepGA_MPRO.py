#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 18:57:14 2021

@author: FCCarvalho
"""
import os
import time
import random as rd
import numpy as np
import pyrosetta
#import modeller
#from modeller.automodel import *
import pymol
from pymol import stored
import joblib
from multiprocessing import Process, Queue
import shutil
import pandas as pd

# Avoiding a very cluttered output
pyrosetta.init("-mute all")
#modeller.log.none()

# Implementation of the four chosen genetic operations
# First operation: Mutation
def mutation(peptide):
    aminoacids = ['A','C','D','E','F','G','H','I','L','M','N','P','Q','R','S','T','V','W','Y'] 
    pep = list(peptide)
    l = len(pep) - 1
    point = rd.randint(0,l)
    while True:
        new_aminoacid = rd.choice(aminoacids)
        if new_aminoacid != pep[point]:
            break
    pep[point] = new_aminoacid
    pep = "".join(pep)
    return pep

# Second operation: Deletion
def deletion(peptide):
    pep = list(peptide)
    l = len(pep) - 1
    point = rd.randint(0,l)
    pep.pop(point)
    pep = "".join(pep)
    return pep

# Third operation: Addition
def insertion(peptide):
    aminoacids = ['A','C','D','E','F','G','H','I','L','M','N','P','Q','R','S','T','V','W','Y']
    pep = list(peptide)
    l = len(pep)
    point = rd.randint(0,l)
    new_aminoacid = rd.choice(aminoacids)
    pep.insert(point, new_aminoacid)
    pep = "".join(pep)
    return pep

# Fourth operation: Crossover
def crossover(peptide1, peptide2, l_limit = 5, u_limit = 30):
    pep1 = list(peptide1)
    pep2 = list(peptide2)
    l1 = len(pep1) - 1
    l2 = len(pep2) - 1
    point1 = rd.randint(1,l1)
    point2 = rd.randint(1,l2)
    if point1 + (l2 - point2 + 1) > u_limit:
        point1 = u_limit - (l2 - point2 + 1)
    elif point1 + (l2 - point2 + 1) < l_limit:
        point1 = l_limit - (l2 - point2 + 1)
    if point2 + (l1 - point1 + 1) > u_limit:
        point2 = u_limit - (l1 - point1 + 1)
    elif point2 + (l1 - point1 + 1) < l_limit:
        point2 = l_limit - (l1 - point1 + 1)
    pep1_l = pep1[:point1]
    pep1_r = pep1[point1:]
    pep2_l = pep2[:point2]
    pep2_r = pep2[point2:]
    pep1 = pep1_l + pep2_r
    pep2 = pep2_l + pep1_r
    pep1 = "".join(pep1)
    pep2 = "".join(pep2)
    return pep1, pep2

def alify(candidate, n = 0):
    
    name = "candidate_" + str(n)
    aliname = name + ".ali"
    code = "cnd" + str(n)
    with open("candidate.ali", 'w') as alifile:
        alifile.write(">P1;cand\n")
        alifile.write("sequence:cand:::::::0.00: 0.00\n")
        alifile.write(candidate)
        alifile.write("*")
        
# Second step: Build the model 
def buildModel(infile = 'candidate', name = 'cand', template ='template', gen = 0):
    # Create an environment
    #env = modeller.Environ()
    
    # Create an alignment
    #aln = modeller.Alignment(env)
    #mdl = modeller.Model(env, file=template)
    aln.append(file='candidate.ali', align_codes='cand')
    aln.append_model(mdl, align_codes='template', atom_files=template +'.pdb')
    aln.align2d()
    aln.write(file='candidate-template.ali', alignment_format='PIR')

    #a = AutoModel(env,
                  #alnfile='candidate-template.ali',      # alignment filename
                  #knowns='template',                      # codes of the templates
                  #sequence='cand')                    # code of the target

    a.very_fast()                 
    a.starting_model = 1
    a.ending_model = 1
    a.make()
    
    # Clear the garbage -> Modeller can't easily personalize file names
    os.remove("cand.D00000001")
    os.remove("cand.ini")
    os.remove("cand.rsr")
    os.remove("cand.sch")
    os.remove("cand.V99990001")
    os.remove("candidate.ali")
    os.remove("candidate-template.ali")
    os.rename("cand.B99990001.pdb", "Results/PDB_predicted/G" + str(gen) + "/" + str(name) + ".pdb")

def calculateFitness_Parallel(candidate, protocol, scorefxn, gen, reps, queue):
    best_pose = None
    best_fitness = 999999999
    scores = list() #Corregir en el principal
    #to_centroid = pyrosetta.SwitchResidueTypeSetMover('centroid')
    for i in range(reps):
        # Load the poses
        comp = pyrosetta.pose_from_pdb("Results/PDB_complex/G" + str(gen) + "/comp(" + candidate + ").pdb")
        #to_centroid.apply(comp)
        # Apply the docking protocol and calculate fitness
        protocol.apply(comp)
        fitness = scorefxn(comp)
        scores.append(fitness)
        if fitness < best_fitness:
            best_pose = comp.clone()
            best_fitness = fitness
    best_pose.dump_pdb("Results/PDB_docking/G" + str(gen) + "/dock(" + candidate + ").pdb")
    queue.put([candidate, scores, best_fitness])

def makeBatch(population, cores, add_fake = False):
    # Population is the list of strings representing the peptides
    # n is the number of peptides per batch
    batches = list()
    x = int(len(population)/cores)
    if len(population) < cores:
        return [population]
    
    for i in range(x):
        if add_fake == True:
            new_batch = ["BJOUXZ"]
        else:
            new_batch = list()
        for j in range(cores):
            pos = i*cores + j
            new_batch.append(population[pos])
        batches.append(new_batch)
        
    if len(population)%cores != 0: # If true, there is an excess of peptides not covered yet
        if add_fake == True:
            final_batch = ["BJOUXZ"]
        else:
            final_batch = list()
        
        for k in range(len(population)%cores):
            pos = x*cores + k
            final_batch.append(population[pos])
        batches.append(final_batch)
    return batches

def makePeptide(candidate, gen=0):
    peptide = pyrosetta.pose_from_sequence(candidate)
    angles = np.random.rand(3,len(candidate))*360 - 180
    for i in range(len(candidate)):
        peptide.set_phi(i+1, angles[0,i])
        peptide.set_psi(i+1, angles[1,i])
        #peptide.set_chi(1, i+1, angles[2,i])   
    src = 'Temp/' + candidate + ".pdb"
    peptide.dump_pdb(src)
    dst = "Results/PDB_predicted/G" + str(gen) + "/" + candidate + ".pdb"
    pymol.cmd.load("template_pos.pdb")
    pymol.cmd.load(src, "new")
    pymol.cmd.align("new", "template_pos")
    pymol.cmd.save(dst, selection = "new")
    pymol.cmd.reinitialize()
    os.remove('Temp/' + candidate + ".pdb")
    makeComplex(candidate, gen)
    
def makeComplex(candidate, gen):
    pymol.cmd.load("target.pdb")
    pymol.cmd.load("Results/PDB_predicted/G" + str(gen) + "/" + candidate + ".pdb", "pep")
    #pymol.cmd.alter('chain A', 'chain="B"')
    pymol.cmd.alter('pep', 'chain="B"')
    #pymol.cmd.alter('chain E', 'chain="A"') # The rename is necessary because Rosetta has a problem in making chain A move
    pymol.cmd.save("Results/PDB_complex/G" + str(gen) + "/comp(" + candidate + ").pdb") # Join both chains in a single PDB as a complex
    pymol.cmd.reinitialize()

def preparePeptide(candidate, gen, Template = "Folded_template"):
    alify(candidate)
    buildModel(name = candidate, template = Template, gen = gen)
    makeComplex(candidate, gen)

def countContacts(candidate_file, target = None, cutoff = 5.0):
    if target == None:
        target = [('20', 'VAL'), ('21', 'THR'), ('23', 'GLY'), ('24', 'THR'), ('25', 'THR'), 
        ('26', 'THR'), ('27', 'LEU'), ('41', 'ALA'), ('42', 'VAL'), ('44', 'CYS'), ('45', 'THR'),
        ('46', 'ALA'), ('49', 'MET'), ('50', 'LEU'), ('52', 'PRO'), ('67', 'LEU'), ('69', 'GLN'), 
        ('118', 'TYR'), ('119', 'ASN'), ('139', 'SER'), ('140', 'PHE'), ('141', 'LEU'), ('142', 'ASN'), 
        ('143', 'GLY'), ('144', 'SER'), ('145', 'CYS'), ('146', 'GLY'), ('163', 'HIS'), ('164', 'HIS'), 
        ('165', 'MET'), ('166', 'GLU'), ('167', 'LEU'), ('168', 'PRO'), ('172', 'HIS'), ('173', 'ALA'), 
        ('186', 'VAL'), ('187', 'ASP'), ('188', 'ARG'), ('189', 'GLN'), ('190', 'THR'), ('191', 'ALA'), 
        ('192', 'GLN'), ('193', 'ALA')]
    pymol.cmd.reinitialize()
    stored.list = []
    pymol.cmd.load(candidate_file)
    pymol.cmd.indicate(''.join(['chain A within ', str(cutoff), ' of chain B']))
    pymol.cmd.iterate('indicate', 'stored.list.append((resi,resn))')
    pymol.cmd.delete('indicate')
    count = 0
    for aa in pd.unique(stored.list):
        if aa in target:
            count += 1
    occupancy = count/len(target)
    return occupancy

def tournament(population, fitness, size = 5, reuse = False):
    # Randomly selects tournament participants among members of the population
    selected = []
    scores = []
    if size >= len(population):
        selected = population.copy()
        scores = fitness.copy()
    else:
        for i in range(size):
            new = rd.choice([p for p in population if p not in selected])
            selected.append(new)
            j = population.index(new)
            scores.append(fitness[j])
    ranking = [x for _,x in sorted(zip(scores,selected))]
    if not reuse:
        j = population.index(ranking[0])
        population.remove(ranking[0])
        fitness.pop(j)
    return ranking[0]

def GA(initial_population, tournament_size = 3, max_generations = 10, pop_size = None, nsims = 20, mutation_rate = 0.20, cross_rate = 0.80, elite = 1, workers=1):
    # Create directories
    if not os.path.isdir("Results/"):
            os.mkdir("Results/")
    if not os.path.isdir("Temp/"):
            os.mkdir("Temp/")
    if not os.path.isdir("Results/PDB_docking/"):
            os.mkdir("Results/PDB_docking/")
    if not os.path.isdir("Results/PDB_predicted/"):
            os.mkdir("Results/PDB_predicted/")
    if not os.path.isdir("Results/PDB_complex/"):
            os.mkdir("Results/PDB_complex/")
            
    removed = []
    # Create results files
    previous_result = False
    initial = 0
    if not os.path.isfile("Results/summary_GA.csv"):
        summary = open("Results/summary_GA.csv", mode="w+") # Create the summary csv file for later visualization of results
        summary.write("generation,best_ind,best_gen_fit,avg_gen_fit,worst_gen_fit,total_time,%surface_occupied,maximum_occupancy" + "\n") # The headers of summary
        summary.close()
    if not os.path.isfile("Results/populations_GA.csv"):
        p_info = open("Results/populations_GA.csv", mode="w+") # Create populations csv file for posterior checks
        p_info.write("generation,genotype,") # The headers of summary
        for i in range(nsims):
            p_info.write("f{},".format(i+1))
        p_info.write("best,occupancy") 
        p_info.close()   
    else:
        previous_result = True
        print("\n"*2, "PREVIOUS RESULTS FILE FOUND. CONTINUING FROM LAST RESULTS.", "\n"*2)
    
    # Calculate population size as the same size of the initial population. This parameter won't change through generations
    if pop_size == None:
        pop_size = len(initial_population)
    else:
        pop_size = 100
    population = initial_population
    
    # Define the protocol
    protocol = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="fa_standard" weights="ref2015.wts"/>
    </SCOREFXNS>
    <MOVERS>
        <FlexPepDock name="ppack"  ppk_only="true"/>
        <FlexPepDock name="fpd" lowres_abinitio="true" pep_refine="true"/>
        <FlexPepDock name="minimize"  min_only="true"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover="ppack"/>
        <Add mover="fpd"/>
        <Add mover="minimize"/>
    </PROTOCOLS>
    <OUTPUT/>
    </ROSETTASCRIPTS>""").get_mover("ParsedProtocol")
    
    #protocol = pyrosetta.rosetta.protocols.docking.DockingLowRes()
    #to_centroid = pyrosetta.SwitchResidueTypeSetMover('centroid')
    
    # Define the scorefxn 
    scorefxn = pyrosetta.create_score_function("ref2015") # This is the standard full_atom scoring from pyrosetta
    #scorefxn = pyrosetta.create_score_function("interchain_cen") # This is the standard centroid scoring from pyrosetta
    
    # Initialize lists to keep records of best, worst and average fitness in each generation
    if previous_result == True:
        pr = pd.read_csv("Results/populations_GA.csv")
        fit_dict = {}
        ocup_dict = {}
        known_genotypes = pr["genotype"].tolist()
        known_scores = pr["best"].tolist()
        known_occupancies = pr["occupancy"].tolist()
        for i in range(len(known_genotypes)):
            fit_dict[known_genotypes[i]]=known_scores[i]
            ocup_dict[known_genotypes[i]]=known_occupancies[i]
        sm = pd.read_csv("Results/summary_GA.csv")
        best_genotypes = sm["best_ind"].tolist()
        best_fitness = sm["best_gen_fit"].tolist()
        average_fitness = sm["avg_gen_fit"].tolist()
        worst_fitness = sm["worst_gen_fit"].tolist()
        initial = sm["generation"].tolist()[-1]
        print("Imported", len(fit_dict), "known peptides")
        print("Starting from generation", initial)
        print("#"*10, "\n"*2)
        max_occupancy = []
        best_occupancy = sm['%surface_occupied'].to_list()
        population = pr.loc[pr['generation'] == initial, ['genotype']]['genotype'].to_list()
    else:
        best_fitness = []
        best_genotypes = []
        worst_fitness = []
        average_fitness = []
        max_occupancy = []
        best_occupancy = []
        fit_dict = {}
        ocup_dict = {}
         
    
    # Calculate fitness for each genotype and fill the records
    for j in range(initial, initial + max_generations):
        print("Beginning generation", j)
        # Prepare the directories
        if not os.path.isdir("Results/PDB_docking/G" + str(j)):
            os.mkdir("Results/PDB_docking/G" + str(j))
        if not os.path.isdir("Results/PDB_predicted/G" + str(j)):
            os.mkdir("Results/PDB_predicted/G" + str(j))
        if not os.path.isdir("Results/PDB_complex/G" + str(j)):
            os.mkdir("Results/PDB_complex/G" + str(j))
            
        # Start timer to check generation duration
        start = time.time() 
        
        # Ensure that the algorithm is not stuck in local minima
        sm = pd.read_csv("Results/summary_GA.csv")
        best_peps = list(sm['best_ind'])
        for b in pd.unique(best_peps):
            if best_peps.count(b) > 3 and pd.unique(sm[sm['best_ind']==b]['best_gen_fit'])[0] < -100:
                removed.append(b)
        
        # Start Queue 
        q = Queue()
        
        # Remove the twins in the generation
        batch_pop = pd.unique(population).tolist()
        
        # Preliminary test to check if there are repeated genotypes from previous generations
        for genotype in reversed(batch_pop):  # If it is resuming simulation, then there will be a duplicate generation
            if genotype in fit_dict.keys() and not previous_result and genotype not in removed:
                # Copy the old pdbs to the new folders and update results
                p_info = open("Results/populations_GA.csv", mode="a") 
                p_info.write("\n"+ str(j) + ","  + genotype)
                for r in range(nsims):
                    p_info.write("," + str(fit_dict[genotype]))
                p_info.write("," + str(fit_dict[genotype]) + "," + str(best_occupancy[-1]))
                p_info.close()
                for r in reversed(range(j)):
                    if os.path.isfile("Results/PDB_docking/G" + str(r) + "/dock(" + genotype + ").pdb"):
                        src1 = "Results/PDB_docking/G" + str(r) + "/dock(" + genotype + ").pdb"
                        dst1 = "Results/PDB_docking/G" + str(j) + "/dock(" + genotype + ").pdb"
                        src2 = "Results/PDB_predicted/G" + str(r) + "/" + genotype + ".pdb"
                        dst2 = "Results/PDB_predicted/G" + str(j) + "/" + genotype + ".pdb"
                        src3 = "Results/PDB_complex/G" + str(r) + "/comp(" + genotype + ").pdb"
                        dst3 = "Results/PDB_complex/G" + str(j) + "/comp(" + genotype + ").pdb"
                        shutil.copyfile(src1, dst1)
                        shutil.copyfile(src2, dst2)
                        shutil.copyfile(src3, dst3)
                        break
                batch_pop.remove(genotype)
                
        # Organize new population in batches
        if previous_result == True:
            batches = []
            occupancies = pr.loc[pr['generation'] == initial, ['occupancy']]['occupancy'].to_list()
        else:
            batches = makeBatch(batch_pop, cores = workers)
            occupancies = list()
        previous_result = False
        # Preparation step is sequential because of Modeller
        z = 0
        
        for batch in batches:
            print("######## BATCH", z, "########\n\n")
            for genotype in batch:
                if genotype != "BJOUXZ": # The fake peptide speeds up simulation since the first process starts alone
                    #if not os.path.isdir("Results/PDB_docking/G0/" + genotype):
                    #    os.mkdir("Results/PDB_docking/G0/" + genotype)
                    #preparePeptide(candidate=genotype, gen=j)
                    makePeptide(candidate=genotype, gen=j)
            
            # Start the process list
            processes = list()
            for genotype in batch: # The docking step is parallel
                p = Process(target=calculateFitness_Parallel, args=(genotype, protocol, scorefxn, j, nsims, q))
                p.start()
                processes.append(p)

            for p in processes: # We can only move forward after all the processes finish
                p.join()
            
            print("##### END OF BATCH", z, "#####")
            print("\n", "#"*10, "Queue Size:", q.qsize(), "#"*10, "\n"*2)
            z += 1
            p_info = open("Results/populations_GA.csv", mode="a")
            for i in range(q.qsize()): # Fetching the scores of each genotype
                [pep, s, bf] = q.get()
                oc = countContacts("Results/PDB_docking/G" + str(j) + "/dock(" + pep + ").pdb")
                occupancies.append(oc)
                scores_txt = [str(score) for score in s]
                scores_txt = ",".join(scores_txt)
                p_info.write("".join(["\n", str(j), ",", pep, ",", scores_txt, ",", str(bf), ",", str(oc)]))
                if oc < 0.3:
                    s = [9999999999, 9999999999] # Make it so the peptide with less than 30% occupancy gets bad scores 
                fit_dict[pep] = min(s) # The best of all the scores is selected as true score
                ocup_dict[pep] = oc
            p_info.close()
            
        fit_list = list()
        for genotype in population:
            fit_list.append(fit_dict[genotype])
            
        if len(batch_pop) > 0: # batch_pop will be an empty list if the simulation is resuming
            best_fitness.append(min(fit_list))
            max_occupancy.append(max(occupancies))
            l = fit_list.index(min(fit_list))
            best_genotypes.append(population[l])
            best_occupancy.append(ocup_dict[best_genotypes[-1]])
            fit_list2 = [f for f in fit_list if f<999999999]
            worst_fitness.append(max(fit_list2))
            average = sum(fit_list2)/len(fit_list2)
            average_fitness.append(average)
        
        # Create next generation
        spots_left = pop_size
        new_pop = []
        if elite > 0:
            spots_left -= elite
            ranking = [x for _,x in sorted(zip(fit_list,population))]
            for n in range(elite):
                best = ranking[n]
                new_pop.append(best)
                # Elite remains able to go through crossover and mutations
        
        # Establish the chances for each operation to happen. Addition, deletion and mutation share mutation chances equally
        operation_weights = [cross_rate, mutation_rate/3, mutation_rate/3, mutation_rate/3]
        while spots_left > 0:
            if spots_left == 1: 
                operation = rd.choice(["insertion", "deletion", "mutation"])
            else:                         
                operation = rd.choices(["crossover","insertion", "deletion", "mutation"], weights=operation_weights)[0]
            if operation == "crossover":
                parent1 = tournament(population,fit_list, tournament_size, reuse = True)
                parent2 = tournament(population,fit_list, tournament_size, reuse = True)
                child1, child2 = crossover(parent1, parent2)
                new_pop.append(child1)
                new_pop.append(child2)
                spots_left -= 2
            elif operation == "insertion":
                chosen = tournament(population,fit_list, tournament_size, reuse = False)
                if len(chosen) < 50:
                    mutant = insertion(chosen)
                else:
                    mutant = deletion(chosen)
                new_pop.append(mutant)
                spots_left -= 1
            elif operation == "deletion":
                chosen = tournament(population,fit_list, tournament_size, reuse = False)
                if len(chosen) > 5:
                    mutant = deletion(chosen)
                else:
                    mutant = insertion(chosen)
                new_pop.append(mutant)
                spots_left -= 1                  
            else:  
                chosen = tournament(population,fit_list, tournament_size, reuse = False)
                mutant = mutation(chosen)
                new_pop.append(mutant)
                spots_left -= 1
        
        end = time.time()
        gen_time = end - start
        if len(batch_pop) > 0: # batch_pop will be an empty list if the simulation is resuming
            summary = open("Results/summary_GA.csv", mode="a") # Open the summary file in append mode
            summary.write(str(j) + "," + str(best_genotypes[-1]) + "," + str(best_fitness[-1]) + "," + str(average) + "," + str(worst_fitness[-1]) + "," + str(gen_time) + "," + str(best_occupancy[-1]) + "," + str(max_occupancy[-1]) + "\n")
            summary.close()
        population = new_pop.copy() #The new population takes over
    return population, best_fitness, worst_fitness, average_fitness, best_genotypes, max_occupancy, best_occupancy

if __name__ == '__main__':
    s = time.time()
    #ranking = pd.read_csv("Propedia_test_Results/populations_GA.csv")
    #initial_population = ranking.sort_values(by = " best")[:90]["genotype"].tolist()
    initial_population = joblib.load('best300.lst')
    p, bf, wf, af, bi, mo, bo = GA(initial_population,max_generations = 3, nsims = 10, pop_size=105, workers=15)
    e = time.time()
    joblib.dump(p, "Results/populations.lst")
    joblib.dump(bf, "Results/best_scores.lst")
    joblib.dump(wf, "Results/worst_scores.lst")
    joblib.dump(af, "Results/avg_scores.lst")
    joblib.dump(bi, "Results/best_genotypes.lst")
    joblib.dump(mo, "Results/max_occupancies.lst")
    joblib.dump(bo, "Results/best_occupancies.lst")
    print('\n\n\n')
    print('FINISHED IN', np.round((e-s)/60,3), 'MINUTES')
