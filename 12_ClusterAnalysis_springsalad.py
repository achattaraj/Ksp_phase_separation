# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 12:26:28 2019

@author: Ani Chattaraj
"""

import re
import os, sys
from decimal import Decimal
import numpy as np
import pandas as pd
from csv import writer
from collections import defaultdict, OrderedDict
from time import time
from glob import glob

#import matplotlib.pyplot as plt
#import multiprocessing as mp


def ProgressBar(jobName, progress, length=40):
    completionIndex = round(progress*length)
    msg = "\r{} : [{}] {}%".format(jobName, "*"*completionIndex + "-"*(length-completionIndex), round(progress*100))
    if progress >= 1: msg += "\r\n"
    sys.stdout.write(msg)
    sys.stdout.flush()

def displayExecutionTime(func):
    """
    This decorator will calculate time needed to execute a function 
    """
    def wrapper(*args, **kwrgs):
        t1 = time()
        func(*args, **kwrgs)
        t2 = time()
        delta = t2 - t1
        if delta < 60:
            print("Execution time : {:.4f} secs".format(delta))
        else:
            t_min, t_sec = int(delta/60), delta%60
            print(f"Execution time : {t_min} mins {t_sec} secs")
    return wrapper


class ReadInputFile:
    
    def __init__(self, txtfile):
        self.txtfile = txtfile
    
    def readFile(self):
        with open(self.txtfile, 'r') as tmpfile:
            txtfile = [line for line in tmpfile.readlines()]
        return txtfile
    
    def getInpath(self):
        mypath = self.txtfile.split("/")[:-1]
        return "/".join(mypath)
    
    def getOutpath(self, statName):
        # statName : Name of the analysis; cluster_stat, e2e_stat etc
        inpath = self.getInpath()
        outpath = inpath + f"/pyStat/{statName}"
        if not os.path.isdir(outpath):
            os.makedirs(outpath)
        return outpath
    
    def getMolecules(self):
        lines = self.readFile()
        molNames, molCounts = [], []
        
        for index, line in enumerate(lines):
            if not re.search('[*]', line):
                if re.search('MOLECULE', line):
                    specs = line.split(':')[-1].split()
                    if len(specs)>2:
                        if int(specs[3]) != 0:  # molecules with non-zero count
                            molNames.append(specs[0].replace('"',''))
                            molCounts.append(int(specs[3]))
                            
        return molNames, molCounts
    
    @staticmethod
    def StringSearch(string, lines):
        value = None
        for line in lines:
            if re.search(string, line):
                value = line.split(':')[-1]
                break
        return value
    
    def getTimeStats(self):
        lines = self.readFile()
        t_tot = self.StringSearch("Total time", lines)
        ts = self.StringSearch("dt", lines)  # ts : time step for simulation
        t_data = self.StringSearch("dt_data", lines)
        t_image = self.StringSearch("dt_image", lines)
        return float(t_tot), float(ts), float(t_data), float(t_image)
    
    def getNumRuns(self):
        numRuns = 0
        lines = self.readFile()
        for line in lines:
            if re.search("\\b" + "Runs" + "\\b", line): # \\b : word boundary; avoids the line 'SimultaneousRuns: 1'
                numRuns = line.split(':')[-1]
                break
        return int(numRuns)



class ClusterAnalysis:

    def __init__(self, simfile, t_total=None, runs=None, withMonomer=False):
        self.simfileObj = ReadInputFile(simfile)

        tf, dt, dt_data, dt_image = self.simfileObj.getTimeStats()
        self.dt = dt_data
        if t_total == None:
            self.t_total = tf
        else:
            self.t_total = t_total
            
        self.withMonomer = withMonomer
        
        timepoints = np.arange(0, self.t_total+self.dt, self.dt)
        self.timePoints = timepoints
        
        numRuns = self.simfileObj.getNumRuns()
        if runs == None:
            self.runs = [run for run in range(numRuns)]
        elif type(runs) is int:
            self.runs = [i for i in range(runs)]
        elif type(runs) is list:
            self.runs = runs
        else:
            print("runs must be an integer or a list of runs")
        

    def __repr__(self):
        simfile = self.simfileObj.txtfile.split('/')[-1]
        info = f"Class : {self.__class__.__name__}\nSystem : {simfile}\nTotal_Time : {self.t_total} seconds\t\tnumRuns : {len(self.runs)}"
        return info

    @staticmethod
    def getSpeciesCount(df, species):
        try:
            # df : pandas dataframe
            col1, col2 = df[0], df[1]
            return [col2[i] for i in range(len(col1)) if col1[i]== species]
        except:
            pass

    def getClusterAverage(self, csvfile):
        # including the monomer
        acs, aco = 0, 0
        clusterList = []
        df = pd.read_csv(csvfile, header=None)
        clusterCount = df[1][0]
        #col1, col2 = df[0], df[1]
        if clusterCount == 0:
            print('No cluster found')
        else:
            oligos = self.getSpeciesCount(df, "Size")
            if len(oligos) == 0:
                acs, aco = 1, 1
            else:
                monomerCount = clusterCount - len(oligos)
                clusterList.extend(oligos)
                clusterList.extend([1]*monomerCount)

                N, TM = len(clusterList), sum(clusterList)
                acs = TM/N # acs : average cluster size
                foTM = {clus: (clusterList.count(clus)*(clus/TM)) for clus in set(clusterList)} #foTM : fraction of total molecules
                aco = sum([cs*f for cs, f in foTM.items()])  # aco : average cluster occupancy

        return acs, aco


    def getOligoAverage(self, csvfile):
        # excluding the monomer
        acs, aco = 0, 0
        df = pd.read_csv(csvfile, header=None)
        oligos = self.getSpeciesCount(df, "Size")
        N, TMC = len(oligos), sum(oligos)  # TMC: total molecules in cluster
        try:
            acs = TMC/N
            foTM = {oligo: (oligos.count(oligo)*(oligo/TMC)) for oligo in set(oligos)}
            aco = sum([cs*f for cs, f in foTM.items()])
        except:
            pass # if there is no cluster > 1, acs = aco = 0
        return acs, aco


    @staticmethod
    def getTimeAverage(ndList, numRuns, numTimePoints):
        ndArray = np.array(ndList).reshape(numRuns, numTimePoints)
        timeAverage = np.mean(ndArray, axis=0, dtype=np.float64) #float64 gives more accurate mean
        return timeAverage

    
    @displayExecutionTime
    def getMeanTrajectory(self):
        
        simObj = self.simfileObj
        inpath = simObj.getInpath()
        if not self.withMonomer:
            outpath = simObj.getOutpath('Cluster_stat_wom')
        else:
            outpath = simObj.getOutpath("Cluster_stat")

        nf = abs(Decimal(str(self.dt)).as_tuple().exponent) # Number of decimal points to format the clusterTime.csv file

        timepoints = self.timePoints
        numRuns = len(self.runs)
        numTimePoints = len(timepoints)
        csvfiles = [inpath + f"/data/Run{run}/Clusters_Time_{tp:.{nf}f}.csv" for run in self.runs for tp in timepoints]
        N = len(csvfiles)
        
        ave_clus_stat = [*map(lambda x: [self.getClusterAverage(x), ProgressBar('Cluster_ave_calc', progress=(csvfiles.index(x)+1)/N)] , csvfiles)]
        ave_stat, _ = list(zip(*ave_clus_stat)) # acs_aco, None
        acs, aco = list(zip(*ave_stat))
        

        print('Done')

        mean_acs =np.mean(np.array(acs).reshape((numRuns, numTimePoints)), axis=0)
        mean_aco =np.mean(np.array(aco).reshape((numRuns, numTimePoints)), axis=0)

        tp_ms = timepoints*1e3 # second to millisecond

        with open(outpath + "/Clustering_dynamics.csv","w", newline='') as outfile:
            wf = writer(outfile,delimiter=',')  # csv writer
            wf.writerow(['Time (ms)','ACS','ACO'])
            wf.writerows(zip(tp_ms, mean_acs, mean_aco))
            outfile.close()



    @staticmethod
    def writeComposition(outpath, compo_dict, molecules):
        d = OrderedDict(sorted(compo_dict.items(), key = lambda x:x[0]))
        with open(outpath + "/Clusters_composition.txt","w") as tmpfile:
            tmpfile.write(f"Cluster Size \t {molecules} : frequency\n\n")
            for k, v in d.items():
                unique_comp = [list(x) for x in set(tuple(x) for x in v)]
                freq = [v.count(comp)/len(v) for comp in unique_comp]
                tmpfile.write(f"  {k}\t\t")
                for cmp, occur in zip(unique_comp, freq):
                    cmp = [str(s) for s in cmp]
                    tmpfile.write(",".join(cmp))
                    tmpfile.write(" : {:.2f}%\t".format(occur*100))
                tmpfile.write("\n\n")

    @staticmethod
    def writeDistribution(outpath, cluster_stat):
        TM_eff, N = sum(cluster_stat), len(cluster_stat)
        unique_clusters = sorted(set(cluster_stat))
        freq = [cluster_stat.count(clus)/N for clus in unique_clusters]
        foTM = [cluster_stat.count(clus)*(clus/TM_eff) for clus in unique_clusters]

        with open(outpath+"/SteadyState_distribution.csv", "w", newline='') as tmpfile:
            wf2 = writer(tmpfile)
            wf2.writerow(['Cluster size','frequency','foTM'])
            wf2.writerows(zip(unique_clusters, freq, foTM))
            tmpfile.close()

    @displayExecutionTime
    def getSteadyStateDistribution(self, SS_timePoints):
        print("Getting steadystate cluster distribution ...")
        simObj = self.simfileObj
        inpath = simObj.getInpath()
        if self.withMonomer:
            outpath = simObj.getOutpath("Cluster_stat")
        else:
            outpath = simObj.getOutpath("Cluster_stat_wom")

        molecules, counts = simObj.getMolecules()
        nf = abs(Decimal(str(self.dt)).as_tuple().exponent) # Number of decimal points to format the clusterTime.csv file
        composition = defaultdict(list) # composition does not track the monomer identity
        #print("nf : ", nf)
        cluster_stat = []
        numRuns = len(self.runs)
        for j, run in enumerate(self.runs):
            for tp in SS_timePoints:
                df = pd.read_csv((inpath + "/data/Run{}/Clusters_Time_{:.{}f}.csv".format(run, tp, nf)), header=None)
                oligos = self.getSpeciesCount(df, "Size")
                mol_in_cluster = [self.getSpeciesCount(df, mol) for mol in molecules]
                for i, clus in enumerate(oligos):
                    comp = [mol_in_cluster[m][i] for m in range(len(mol_in_cluster))]
                    composition[clus].append(comp)

                if not self.withMonomer:
                    cluster_stat.extend(oligos)
                else:
                    monomerCount = df[1][0] - len(oligos)
                    monomers = [1]*monomerCount
                    cluster_stat.extend(monomers)
                    cluster_stat.extend(oligos)

            ProgressBar("Progress", (j+1)/numRuns)

        self.writeComposition(outpath, composition, molecules)
        self.writeDistribution(outpath, cluster_stat)

        with open(outpath + "/Sampling_stat.txt","w") as tmpfile:
            ss_tp1000 = [(t*1e3) for t in SS_timePoints]
            tmpfile.write("Number of runs : {}\n\n".format(len(self.runs)))
            tmpfile.write(f"Steady state time points (ms): {ss_tp1000}")



if __name__ == "__main__":
    
    file = "Z:/01_PROJECTS/00_PS_principles/01_ssalad_manuscript_models/00_Reference_system/4v_4v_flex_9nm_linker_SIMULATIONS/4v_flex_config03_CB_25uM_SIM_FOLDER/4v_flex_config03_CB_25uM_SIM.txt"
    runs = [1,2,3]
    ca = ClusterAnalysis(file, t_total=None, runs=runs, withMonomer=True)
    print(ca)
    ca.getMeanTrajectory()
    ca.getSteadyStateDistribution(SS_timePoints=[0.05])
    
    


    
    
