# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 17:33:23 2020

@author: Ani Chattaraj
"""

import numpy as np
from glob import glob
import sys, os
import re
from collections import defaultdict, OrderedDict
import csv


class NFSim_output_analyzer:
    def __init__(self, path):
        '''
        Parameters
        ----------
        path : File String
            DESCRIPTION: location of the source directory containing gdat files

        Returns
        -------
        None.

        '''
        self.path = path
    def __repr__(self):
        simfile = self.path.split('/')[-1]
        gfiles = glob(self.path + "\*.gdat")
        #sfiles = glob(self.path + "\*.species")
        info = f"\n***** // ***** \nClass : {self.__class__.__name__}\nSystem : {simfile}\nTotal Trajectories : {len(gfiles)}\n"
        return info

    #@displayExecutionTime
    def process_gdatfiles(self):
        '''

        Computes Mean observable counts over multiple trajectories

        '''
        gfiles = glob(self.path + "\*.gdat")
        if len(gfiles) == 0:
            print('No gdat files found; quitting calculation ...')
            sys.exit()

        '''I use a test gdat file to extract the array dimension
            and name of the observables used in the model'''

        test_gf = gfiles[0]
        N_tp, N_obs = np.loadtxt(test_gf).shape # number of timepoints and observables

        with open(test_gf,'r') as tmpf:
            obs_names = tmpf.readline().split()[1:]
            obs_names = '\t'.join(obs_names)


        '''The temporary matrix would store the data from multiple trajectories
            and perform the average'''

        tmp_matrix = np.empty(shape=(len(gfiles),N_tp, N_obs), dtype=float)
        N_gf = len(gfiles)

        for i, gf in enumerate(gfiles):
            data = np.loadtxt(gf)
            tmp_matrix[i] = data
            self.ProgressBar('Processing gdat_files', (i+1)/N_gf)

        mean_obs = np.mean(tmp_matrix, axis=0)
        outpath = self.getOutpath()
        np.savetxt(outpath + "\Mean_Observable_Counts.txt", mean_obs, header=obs_names, fmt='%.6e')

    #@displayExecutionTime
    def process_speciesfiles(self, molecules=[]):

        '''

        molecules = List of molecules used in the model
        -------

        Computes distribution of molecular clusters and their compositions

        '''
        sfiles = glob(self.path + "/*.species")
        flatten_ = lambda myList: [item for sublist in myList for item in sublist]
        cs_stat, comp_stat = defaultdict(list), defaultdict(list)
        N_sp = len(sfiles)

        for i, sf in enumerate(sfiles):
            cs, comp = self.collect_clusters(sf, molecules)
            #print('cs = ', cs)
            for size, count in cs.items():
                cs_stat[size].append(count)
            for size, composition in comp.items():
                comp_stat[size].append(composition)
            self.ProgressBar('Processing species_files', (i+1)/N_sp)

        cs_stat = {k: sum(v) for k, v in cs_stat.items()}
        comp_stat = {k: flatten_(v) for k, v in comp_stat.items()}

        outpath = self.getOutpath()

        self.writeComposition(outpath, comp_stat, molecules)
        self.writeDistribution(outpath, cs_stat)

    def getOutpath(self):
        outpath = self.path + "/pyStat"
        if not os.path.isdir(outpath):
            os.makedirs(outpath)
        return outpath

    @staticmethod
    def collect_clusters(speciesFile, molecules):
        '''
        Parameters
        ----------
        speciesFile : File String
            DESCRIPTION: Speciesfile containing all the molecular species
        molecules : List of String
            DESCRIPTION: List of molecules used in the model

        Returns
        -------
        A pair of defaultdicts; one with cluster size distribution
        and another with the corresponding compositions of the clusters

        '''
        try:
            with open(speciesFile, 'r') as tf:
                currentFrame = tf.readlines()[2:] # to avoid first two warning lines
        except:
            print("File missing: ", speciesFile)
            sys.exit()
        else:
            clus_stat = defaultdict(list)
            comp_stat = defaultdict(list)
            for line in currentFrame:
                if not (line == '\n' or re.search('Time', line) or re.search('Sink', line) or re.search('Source', line)):
                    cluster, count = line.split()
                    comp = tuple([cluster.count(mol) for mol in molecules])
                    cs = len(cluster.split('.'))
                    if cs == 0:
                        cs == 1  # monomer does not have bonds (.)

                    clus_stat[cs].append(int(count))
                    comp_stat[cs].append(comp)
        clus_stat = {k: sum(v) for k,v in clus_stat.items()}
        return clus_stat, comp_stat


    @staticmethod
    def writeDistribution(outpath, cluster_stat):
        '''
        Parameters
        ----------
        outpath : File String
            DESCRIPTION: Location of the output files
        cluster_stat : Defaultdict
            DESCRIPTION: Dictionary with {keys, values} = {cluster size, occurence}

        Returns
        -------
        None.
        '''
        cluster_stat = OrderedDict(sorted(cluster_stat.items(), key = lambda x:x[0]))
        TC = sum(cluster_stat.values()) # total counts
        TM = sum([k*v for k,v in cluster_stat.items()])  # total molecules
        #print('TM = ', TM, ' TC = ',  TC)
        foTM = {cs: count*(cs/TM) for cs,count in cluster_stat.items()} # fraction of total molecules
        occurence = {cs: count/TC for cs, count in cluster_stat.items()}

        with open(outpath + "/Cluster_frequency.csv","w", newline='') as tmpfile:
            wf = csv.writer(tmpfile)
            wf.writerow(['Cluster size','counts'])
            wf.writerows(zip(cluster_stat.keys(),cluster_stat.values()))

        with open(outpath+"/SteadyState_distribution.csv", "w", newline='') as tmpfile:
            wf2 = csv.writer(tmpfile)
            wf2.writerow(['Cluster size','frequency','foTM'])
            wf2.writerows(zip(cluster_stat.keys(), occurence.values(), foTM.values()))

    @staticmethod
    def writeComposition(outpath, compo_dict, molecules):
        '''
        Parameters
        ----------
        outpath : File String
            DESCRIPTION: Location of the output files
        compo_dict : Defaultdict
            DESCRIPTION: Dictionary with {keys, values} = {cluster size, compositions}
        molecules : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        d = OrderedDict(sorted(compo_dict.items(), key = lambda x:x[0]))
        with open(outpath + "/Clusters_composition.txt","w") as tmpfile:
            tmpfile.write(f"Cluster Size \t {molecules} : frequency\n\n")

            for k, v in d.items():
                unique_comp = set(v)
                freq = [v.count(uc)/len(v) for uc in unique_comp]
                tmpfile.write(f"  {k}\t\t")
                for cmp, occur in zip(unique_comp, freq):
                    cmp = [str(s) for s in cmp]
                    tmpfile.write(",".join(cmp))
                    tmpfile.write(" : {:.2f}%\t".format(occur*100))

                tmpfile.write("\n\n")

    @staticmethod
    def ProgressBar(jobName, progress, length=40):

        '''
        Parameters
        ----------
        jobName : string
            Name of the job given by user.

        progress : float
            progress of the job to be printed as percentage.

        length : interger
            prints the length of the progressbar. The default is 40.

        Returns
        -------
        None.
        '''
        completionIndex = round(progress*length)
        msg = "\r{} : [{}] {}%".format(jobName, "*"*completionIndex + "-"*(length-completionIndex), round(progress*100))
        if progress >= 1: msg += "\r\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    @staticmethod
    def displayExecutionTime(func):
        """
        This decorator (function) will calculate the time needed to execute a task
        """
        def wrapper(*args, **kwrgs):
            t1 = time.time()
            func(*args, **kwrgs)
            t2 = time.time()
            delta = t2 - t1
            if delta < 60:
                print("Execution time : {:.4f} secs".format(delta))
            else:
                t_min, t_sec = int(delta/60), delta%60
                print(f"Execution time : {t_min} mins {t_sec} secs")

        return wrapper


if __name__ == '__main__':
    
    mypath = "Z:/NFSim/4v_4v_FTC_60uM_sample"

    molecules = ['poly_A', 'poly_B']
    #molecules = ['Nephrin', 'Nck', 'NWASP']
    
    nfs_obj = NFSim_output_analyzer(mypath)
    print(nfs_obj)
    nfs_obj.process_gdatfiles()
    nfs_obj.process_speciesfiles(molecules)


   

       
