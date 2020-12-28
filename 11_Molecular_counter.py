# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 17:25:41 2020

@author: Ani Chattaraj
"""

import numpy as np
import os, re
from progress_bar import updateCalc
from time import time
from csv import writer, reader
from AllBondPlot import PlotAllBond
import matplotlib.pyplot as plt
import pandas as pd


class MoleculeCounter:
    
    def __init__ (self, inpath, numRuns):
        self.inpath = inpath
        if type(numRuns) == list:
            self.runs = numRuns
        if type(numRuns) == int:
            self.runs = [_ for _ in range(numRuns)]
    
    
    def getOutpath(self):
        outpath = self.inpath + "/pyStat/Bond_Stat"
        if not os.path.isdir(outpath):
            os.makedirs(outpath)
        return outpath
    
    def getBoxSize(self):
        txtfile = self.inpath.split('/')[-1].replace('_FOLDER','.txt')
        with open(self.inpath + "/" + txtfile, 'r') as tf:
            for line in tf.readlines():
                if re.search('L_x', line):
                    Lx = float(line.split(':')[-1].strip())
                if re.search('L_y', line):
                    Ly = float(line.split(':')[-1].strip())
                if re.search('L_z_in', line):
                    Lz = float(line.split(':')[-1].strip())
        return (Lx*1e3,Ly*1e3,Lz*1e3)  # in nm
    
    @staticmethod
    def getMoleculeCount(file, mols):
        mol_stat = []
        df = pd.read_csv(file, header=0)
        for mol in mols:
            mol_stat.append(list(df[mol]))
        
        return np.asarray(mol_stat)
    
    def getTimePoints(self):
        file = None
        for run in self.runs:
            testfile = self.inpath + f"/data/Run{run}/FullBondData.csv"
            if os.path.isfile(testfile):
                file = testfile
                break
            else:
                pass
        else:
            print('No run is complete yet')
        df = pd.read_csv(file)
       
        return df['Time']
    
    def getMoleculeStat(self, molecules):
        molCount = []
        
        tp = self.getTimePoints()
        
        lx,ly,lz = self.getBoxSize()
        sysVol = lx*ly*lz
        
        for run in self.runs:
            try:
                file = self.inpath + "/data/Run{}/FullCountData.csv".format(run)
                mols = self.getMoleculeCount(file, molecules)
                molCount.append(mols)
                
            except:
                print(f'Run{run} is still not complete!' )
                pass
    
        molCount_arr = np.asarray(molCount)
        
        _, numMols, _ = molCount_arr.shape # shape : numRuns, numMolecules, numTimepoints
        
        m_mols = np.mean(molCount_arr, axis=0) # axis = 0 gives average molecular counts over multiple trajectories
        
        factor = 6.023 * 1e-7  # converts uM to molecules/nm3
        
        df = pd.DataFrame({"Time":[_ for _ in tp]})
        
        for i in range(numMols):
            df[molecules[i] + ' (uM)'] = [_ for _ in m_mols[i]/(factor*sysVol)]
            
        df.to_csv(self.getOutpath() + "/Average_Molecular_Concentration.csv", sep=',')
        


if __name__ == '__main__':
    
    path = 'Z:/01_PROJECTS/00_PS_principles/01_ssalad_manuscript_models/00_Reference_system/4v_4v_flex_9nm_linker_SIMULATIONS/4v_flex_config03_CB_25uM_SIM_FOLDER'
    
    #mols = ['FREE poly_SH3', 'FREE poly_PRM']
    mols = [' TOTAL poly_A', 'TOTAL poly_B']
    
    Runs = 100
    
    mc = MoleculeCounter(path, Runs)
    mc.getMoleculeStat(mols)
    
            
            