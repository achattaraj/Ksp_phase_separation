# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 17:47:25 2020

@author: Ani Chattaraj
"""
from collections import namedtuple

'''
This object is used to specify molecular components.
name - name of the molecule
sites - number (name) of sites 
in_rate - creation rate of that molecule (unit: molecules/s) 
out_rate - decay rate of that molecule (unit: 1/s)
'''
molecule = namedtuple('molecule',['name','sites','in_rate','out_rate'])


class ModelBuilder:

    def __init__(self, model_name, moleculeObjs, outpath,
                 initialCounts=[0,0], kd=1, koff=1,
                 timepoints=[], utl=2):
        self.modelName = model_name
        self.outpath = outpath
        self.molObj = moleculeObjs  # Reactive molecules with number of sites and creation and decay rates
        self.IC = initialCounts     # Initial number of molecules
        self.kd = kd        # Dissociation constant for Individual binding (molecules)
        self.koff = koff    # off rate constant (1/s)
        self.tp = timepoints # timepoints = [t_end, numSteps]
        self.utl = utl      # universal traverse limit

    @staticmethod
    def format_molecule(molObj):
        molString = [molObj.name,'(']
        L = len(molObj.sites)
        for i,site in enumerate(molObj.sites):
            if i != L-1:
                molString.append(site)
                molString.append(',')
            else:
                molString.append(site)
        molString.append(')')
        return ''.join(molString)

    def writeModel(self):
        modelString = f'{self.outpath}/{self.modelName}.bngl'
        mf = open(modelString, 'w')
        print('writing model... \n', modelString)
        mf.write("begin model\n\n")

        mf.write('begin parameters\n')
        mf.write(f'kd\t{self.kd}\n')
        kon=self.koff/self.kd
        mf.write(f'kon\t{kon:.6f}\n')
        mf.write(f'koff\t{self.koff}\n')
        mf.write('end parameters\n\n')

        mf.write('begin molecule types\n')
        mf.write ('Source(s)\n')
        #mf.write ('Sink(s)\n')
        for mol in self.molObj:
            mol_site = self.format_molecule(mol)
            mf.write(f'{mol_site}\n')
        mf.write ('end molecule types\n\n')


        mf.write('begin seed species\n')
        mf.write('1 Source(s) 1.0\n')
        #mf.write('2 Sink(s) 0.0\n')
        speciesCount = 3
        for mol, ic in zip(self.molObj, self.IC):
            mol_site = self.format_molecule(mol)
            mf.write(f'{speciesCount} {mol_site} {ic}\n')
            speciesCount += 1
        mf.write('end seed species\n\n')

        molA, molB = self.molObj[0], self.molObj[1]

        mf.write('begin observables\n')
        for mol in self.molObj:
            free_mol = self.format_molecule(mol)
            mf.write(f'Molecules tot_{mol.name} {mol.name}()\n')
            mf.write(f'Molecules free_{mol.name} {free_mol}\n')

        mf.write('end observables\n\n')



        mf.write('begin reaction rules\n')
        for mol in self.molObj:
            mol_site = self.format_molecule(mol)
            mf.write(f'create_{mol.name}: Source(s) -> Source(s) + {mol_site}\t{mol.in_rate}\n')
            mf.write(f'decay_{mol.name}:  {mol_site} -> 0\t{mol.out_rate}\n')
        rxn_id = 1
        for sa in molA.sites:
            for sb in molB.sites:
                mf.write(f'Rule_{rxn_id}:  {molA.name}({sa}) + {molB.name}({sb}) <-> {molA.name}({sa}!1).{molB.name}({sb}!1)  kon, koff\n')
                rxn_id += 1

        mf.write('end reaction rules\n\n')
        mf.write('end model\n\n')

        brak_open, brak_close = '{', '}'
        mf.write(f'simulate_nf({brak_open}t_end=>{self.tp[0]},n_steps=>{self.tp[1]}, param=>"-v -utl {self.utl}"{brak_close})\n')

        mf.close()
        print('...done writing!')

'''
Model creation for NFSim
'''
outpath = "C:\\Users\\achattaraj\\Dropbox\\2019-05-Visual\\BioNetGen-2.4.0-Win\\PS_principles"

model_name = '4v_4v_FTC_60uM'

nck_sites = [f'a{i}' for i in range(1,5)]
nw_sites = [f'b{i}' for i in range(1,5)]
#print(tuple(sites_1))
nck = molecule(name='poly_A', sites=nck_sites, in_rate=0, out_rate=0)
nwasp = molecule(name='poly_B', sites=nw_sites, in_rate=0, out_rate=0)

mols = [nck, nwasp]

initialCounts = [600, 600]

utl = 2 # universal traverse limit

mb = ModelBuilder(model_name, mols, outpath, initialCounts,
                  kd=3500, koff=100, timepoints=[0.5,100], utl=utl)
mb.writeModel()





