import pmx
from pmx.utils import create_folder
from pmx import *
import sys
import os,shutil
import re
import subprocess
import glob
import random
import pandas as pd
import numpy as np

class AAtutorial: # name modified
    """Class contains parameters for setting up free energy calculations

    Parameters
    ----------
    ...

    Attributes
    ----------
    ....

    """

    def __init__(self, **kwargs):
        
        # set gmxlib path
        gmx.set_gmxlib()
        
        # the results are summarized in a pandas framework
        self.resultsAll = pd.DataFrame()
        self.resultsSummary = pd.DataFrame()
        
        # paths
        self.workPath = './'
        self.pdbfile = 'inputs/protein.pdb'
        self.mdpPath = '{0}/mdp'.format(self.workPath)

        self.edges = {}
        self.resids = {}
        
        # parameters for the general setup
        self.replicas = 5
        self.simTypes = ['em','eq', 'transitions'] # modified
        self.states = ['stateA','stateB']
                
        # simulation setup
        self.ff = 'amber99sb-star-ildn-mut'
        self.boxshape = 'dodecahedron'
        self.boxd = 1.0
        self.water = 'tip3p'
        self.conc = 0.15
        self.pname = 'NaJ'
        self.nname = 'ClJ'
        
        # job submission params
        self.JOBqueue = 'SGE' # could be SLURM
        self.JOBsimtime = 24 # hours
        self.JOBnotr = 100
        self.JOBsimcpu = 8 # CPU default
        self.JOBbGPU = True
        self.JOBmodules = []
        self.JOBsource = []
        self.JOBexport = []
        self.JOBgmx = 'gmx mdrun'
        self.JOBpartition = ''

        for key, val in kwargs.items():
            setattr(self,key,val)
            
    def prepareFreeEnergyDir( self ):
        
        
        # read edges (directly or from a file)
        self._read_edges()
        
        # read mdpPath
        self.mdpPath = self._read_path( self.mdpPath )
        
        # workpath
        self.workPath = self._read_path( self.workPath )
        create_folder( self.workPath )
        
        # create folder structure
        self._create_folder_structure( )
        
        # print summary
        self._print_summary( )

        # print folder structure
        self._print_folder_structure( )
        
        print('DONE')
        
        
    # _functions to quickly get a path at different levels, e.g wppath, edgepath... like in _create_folder_structure
    def _get_specific_path( self, edge=None, state=None, r=None, sim=None ):
        if edge==None:
            return(self.workPath)       
        edgepath = '{0}/{1}'.format(self.workPath,edge)
        
        if state==None:
            return(edgepath)
        statepath = '{0}/{1}'.format(edgepath,state)
        
        if r==None:
            return(statepath)
        runpath = '{0}/run{1}'.format(statepath,r)
        
        if sim==None:
            return(runpath)
        simpath = '{0}/{1}'.format(runpath,sim)
        return(simpath)
                
    def _read_path( self, path ):
        return(os.path.abspath(path))
                
    def _read_edges( self ):
        # read from file
        try:
            if os.path.isfile( self.edges ):
                self._read_edges_from_file( self )
        # edge provided as an array
        except: 
            foo = {}
            for re,e in enumerate(self.edges):
                key = '{0}{1}{2}'.format(e[0],self.resids[re],e[1])
                foo[key] = e
            self.edges = foo
            
    def _read_edges_from_file( self ):
        self.edges = 'Edges read from file'
        
        
    def _create_folder_structure( self, edges=None ):
        # edge
        if edges==None:
            edges = self.edges        
        for edge in edges:
            print(edge)            
            edgepath = '{0}/{1}'.format(self.workPath,edge)
            create_folder(edgepath)

            # stateA/stateB
            for state in self.states:
                statepath = '{0}/{1}'.format(edgepath,state)
                create_folder(statepath)

                # run1/run2/run3
                for r in range(1,self.replicas+1):
                    runpath = '{0}/run{1}'.format(statepath,r)
                    create_folder(runpath)

                    # em/eq_posre/eq/transitions
                    for sim in self.simTypes:
                        simpath = '{0}/{1}'.format(runpath,sim)
                        create_folder(simpath)

    def _print_summary( self ):
        print('\n---------------------\nSummary of the setup:\n---------------------\n')
        print('   workpath: {0}'.format(self.workPath))
        print('   pdb file to use: {0}'.format(self.pdbfile))
        print('   mdp path: {0}'.format(self.mdpPath))
        print('   number of replicase: {0}'.format(self.replicas))
        print('   mutations:')
        for e in self.edges.keys():
            print('        {0}'.format(e))

    def _print_folder_structure( self ):
        print('\n---------------------\nDirectory structure:\n---------------------\n')
        print('{0}/'.format(self.workPath))
        print('|')
        print('|--X2Y')
        print('|--|--stateA')
        print('|--|--|--run1/2/3')
        print('|--|--|--|--em/eq/transitions')
        print('|--|--stateB')
        print('|--|--|--run1/2/3')
        print('|--|--|--|--em/eq/transitions')
        print('|--X2Y')

    def _be_verbose( self, process, bVerbose=False ):
        out = process.communicate()            
        if bVerbose==True:
            printout = out[0].splitlines()
            for o in printout:
                print(o)
        # error is printed every time                  
        printerr = out[1].splitlines()                
        for e in printerr:
            print(e)              

            
    def hybrid_structure_topology( self, edges=None, bVerbose=False ):
        print('----------------------------------')
        print('Creating hybrid structure/topology')
        print('----------------------------------')

        if edges==None:
            edges = self.edges        
        for i, edge in enumerate(edges):
            print(edge)
            outpath = self._get_specific_path(edge=edge)
            protpdb = '{0}'.format(self.pdbfile)
            #
            m = Model(protpdb, rename_atoms=True)
            m2 = mutate(m=m, mut_resid=self.resids[i], mut_resname=self.edges[edge][1], ff=self.ff)
            m2.write('{0}/mutant.pdb'.format(outpath))
            gmx.pdb2gmx(f='{0}/mutant.pdb'.format(outpath), o='{0}/conf.pdb'.format(outpath), p='{0}/topol_prev.top'.format(outpath), ff=self.ff, water='tip3p', other_flags=' -i {0}/posre.itp'.format(outpath))
            topol = Topology('{0}/topol_prev.top'.format(outpath), ff=self.ff)
            pmxtop, _ = gen_hybrid_top(topol)
            pmxtop.write('{0}/topol.top'.format(outpath), scale_mass=0.33)

        print('DONE')
            
            
            

        
    def _clean_backup_files( self, path ):
        toclean = glob.glob('{0}/*#'.format(path)) 
        for clean in toclean:
            os.remove(clean)        
    
    def boxWaterIons( self, edges=None, bBoxProt=True, bWatProt=True, bIonProt=True):
        print('----------------')
        print('Box, water, ions')
        print('----------------')
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            print(edge)            
            outPath = self._get_specific_path(edge=edge)

            # box protein
            if bBoxProt==True:            
                inStr = '{0}/conf.pdb'.format(outPath)
                outStr = '{0}/box.pdb'.format(outPath)
                gmx.editconf(inStr, o=outStr, bt=self.boxshape, d=self.boxd, other_flags='')

            # water protein
            if bWatProt==True:            
                inStr = '{0}/box.pdb'.format(outPath)
                outStr = '{0}/water.pdb'.format(outPath)
                top = '{0}/topol.top'.format(outPath)
                gmx.solvate(inStr, cs='spc216.gro', p=top, o=outStr)                
            
            # ions protein
            if bIonProt:
                inStr = '{0}/water.pdb'.format(outPath)
                outStr = '{0}/ions.pdb'.format(outPath)
                mdp = '{0}/em_l0.mdp'.format(self.mdpPath)
                tpr = '{0}/tpr.tpr'.format(outPath)
                top = '{0}/topol.top'.format(outPath)
                mdout = '{0}/mdout.mdp'.format(outPath)
                gmx.grompp(f=mdp, c=inStr, p=top, o=tpr, maxwarn=4, other_flags=' -po {0}'.format(mdout))
                gmx.genion(s=tpr, p=top, o=outStr, conc=self.conc, neutral=True, 
                      other_flags=' -pname {0} -nname {1}'.format(self.pname, self.nname))                
                     
            # clean backed files
            self._clean_backup_files( outPath )
        print('DONE')
            
    def _prepare_single_tpr( self, simpath, toppath, state, simType, empath=None, frameNum=0 ):
        
        mdpPrefix = ''
        if simType=='em':
            mdpPrefix = 'em'
        elif simType=='eq':
            mdpPrefix = 'eq'
        elif simType=='transitions':
            mdpPrefix = 'ti'        
            
        top = '{0}/topol.top'.format(toppath)
        tpr = '{0}/tpr.tpr'.format(simpath)
        mdout = '{0}/mdout.mdp'.format(simpath)
        # mdp
        if state=='stateA':
            mdp = '{0}/{1}_l0.mdp'.format(self.mdpPath,mdpPrefix)
        else:
            mdp = '{0}/{1}_l1.mdp'.format(self.mdpPath,mdpPrefix)
        # str
        if simType=='em':
            inStr = '{0}/ions.pdb'.format(toppath)
        elif simType=='eq':
            inStr = '{0}/confout.gro'.format(empath)
        elif simType=='transitions':
            inStr = '{0}/frame{1}.gro'.format(simpath,frameNum)
            tpr = '{0}/ti{1}.tpr'.format(simpath,frameNum)
            
        gmx.grompp(f=mdp, c=inStr, p=top, o=tpr, maxwarn=4, other_flags=' -po {0}'.format(mdout))
        self._clean_backup_files( simpath )
                    
         
    def prepare_simulation( self, edges=None, simType='em'):
        print('-----------------------------------------')
        print('Preparing simulation: {0}'.format(simType))
        print('-----------------------------------------')
        
        mdpPrefix = ''
        if simType=='em':
            mdpPrefix = 'em'
        elif simType=='eq':
            mdpPrefix = 'eq'
        elif simType=='transitions':
            mdpPrefix = 'ti'
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            print(edge)
            toppath = self._get_specific_path(edge=edge)
            
            for state in self.states:
                for r in range(1,self.replicas+1):
                    simpath = self._get_specific_path(edge=edge,state=state,r=r,sim=simType)
                    empath = self._get_specific_path(edge=edge,state=state,r=r,sim='em')
                    self._prepare_single_tpr( simpath, toppath, state, simType, empath )
        print('DONE')
        

    def _run_mdrun( self, tpr=None, ener=None, confout=None, mdlog=None, 
                    cpo=None, trr=None, xtc=None, dhdl=None, bVerbose=False):
        # EM
        if xtc==None:
            process = subprocess.Popen(['gmx','mdrun',
                                '-s',tpr,
                                '-e',ener,
                                '-c',confout,
                                '-o',trr,                                        
                                '-g',mdlog],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
            self._be_verbose( process, bVerbose=bVerbose )                    
            process.wait()           
        # other FE runs
        else:
            process = subprocess.Popen(['gmx','mdrun',
                                '-s',tpr,
                                '-e',ener,
                                '-c',confout,
                                '-dhdl',dhdl,
                                '-x',xtc,
                                '-o',trr,
                                '-cpo',cpo,                                        
                                '-g',mdlog],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
            self._be_verbose( process, bVerbose=bVerbose )                    
            process.wait()           
            

    def run_simulation_locally( self, edges=None, simType='em', bVerbose=False ):
        print('-------------------------------------------')
        print('Run simulation locally: {0}'.format(simType))
        print('-------------------------------------------')
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            
            for state in self.states:
                for r in range(1,self.replicas+1):            
                    print('Running: PROT {0} {1} run{2}'.format(edge,state,r))
                    simpath = self._get_specific_path(edge=edge,state=state,r=r,sim=simType)
                    tpr = '{0}/tpr.tpr'.format(simpath)
                    ener = '{0}/ener.edr'.format(simpath)
                    confout = '{0}/confout.gro'.format(simpath)
                    mdlog = '{0}/md.log'.format(simpath)
                    trr = '{0}/traj.trr'.format(simpath)
                    self._run_mdrun(tpr=tpr,trr=trr,ener=ener,confout=confout,mdlog=mdlog,bVerbose=bVerbose)
                    self._clean_backup_files( simpath )
        print('DONE')

    def _create_jobscript(self,fname,jobname,modules,gmx,partition,simpath,indx):
        fp = open(fname,'w')
        fp.write('#!/bin/bash\n')
        fp.write('#SBATCH --job-name={}\n'.format(jobname))
        fp.write('#SBATCH --get-user-env\n')
        fp.write('#SBATCH --nodes=1\n')
        fp.write('#SBATCH --ntasks-per-node={}\n'.format(indx[1]-indx[0]))
        fp.write('#SBATCH --cpus-per-task={}\n'.format(self.JOBntomp))
        fp.write('#SBATCH --mail-type=none\n')
        fp.write('#SBATCH --time=24:00:00\n')
        fp.write('#SBATCH --partition={}\n'.format(partition))
        if self.JOBbGPU == True:
            fp.write('#SBATCH --gres=gpu:1\n')
        fp.write('#SBATCH --exclusive\n')
        fp.write('\n')

        for source in self.JOBsource:
            fp.write('source {}\n'.format(source))
        fp.write('\n')
        for module in modules:
            fp.write('module load {}\n'.format(module))
        fp.write('export GMXRUN="{}"\n'.format(gmx))
        fp.write('\n')

        fp.write('\n')
        fp.write('cd $TMPDIR\n')
        fp.write('for i in {{{}..{}}}; do \n'.format(indx[0]+1, indx[1]))
        fp.write('mkdir ti$i\n')
        fp.write('cp {}/ti$i.tpr ti$i/tpr.tpr\n'.format(simpath))
        fp.write('done\n')
        fp.write('$GMXRUN -s tpr.tpr -multidir ti{{{}..{}}}\n'.format(indx[0]+1, indx[1]))
        fp.write('\n')
        fp.write('for i in {{{}..{}}}; do \n'.format(indx[0]+1, indx[1]))
        fp.write('cp ti$i/dhdl.xvg {}/dhdl$i.xvg\n'.format(simpath))
        fp.write('done\n')
        fp.close()


        fp.close()
    def prepare_jobscripts( self, edges=None, simType='em'):
        print('---------------------------------------------')
        print('Preparing jobscripts for: {0}'.format(simType))
        print('---------------------------------------------')
        
        jobfolder = '{0}/{1}_jobscripts'.format(self.workPath,simType)
        os.system('mkdir {0}'.format(jobfolder))
        
        if edges==None:
            edges = self.edges
            
        counter = 0
        for edge in edges:
            
            for state in self.states:
                for r in range(1,self.replicas+1):            
                    simpath = self._get_specific_path(edge=edge,state=state,r=r,sim=simType)
                    jobfile = '{0}/jobscript{1}'.format(jobfolder,counter)
                    jobname = 'j{0}_{1}_{2}_{3}'.format(edge,state,r,simType)
                    job = pmx.jobscript.Jobscript(fname=jobfile,
                                    queue=self.JOBqueue,simcpu=self.JOBsimcpu,
                                    jobname=jobname,modules=self.JOBmodules,source=self.JOBsource,
                                    gmx=self.JOBgmx,partition=self.JOBpartition,bGPU=self.JOBbGPU)
                    cmd1 = 'cd {0}'.format(simpath)
                    cmd2 = '$GMXRUN -s tpr.tpr'
                    job.cmds = [cmd1,cmd2]
                    if simType=='transitions':
                        if self.JOBnotr%self.JOBsimcpu==0:
                            njobs = int(self.JOBnotr/self.JOBsimcpu)
                        else:
                            njobs = int(self.JOBnotr/self.JOBsimcpu)+1
                        for nj in range(njobs):
                            if nj==njobs-1:
                                indx=[int(nj*self.JOBsimcpu), self.JOBnotr]
                            else:
                                indx = [int(nj*self.JOBsimcpu), int((nj+1)*self.JOBsimcpu)]

                            jobfile = '{0}/jobscript{1}'.format(jobfolder,counter)
                            jobname = 'prot_{0}'.format(counter)
                            self._create_jobscript(fname=jobfile,jobname=jobname,modules=self.JOBmodules,gmx=self.JOBgmx,partition=self.JOBpartition,simpath=simpath,indx=indx)
                            counter+=1
                    else:
                        job.create_jobscript()
                        counter+=1
                        
        #######
        self._submission_script( jobfolder, counter, simType )
        print('DONE')
        
    def _submission_script( self, jobfolder, counter, simType='eq' ):
        fname = '{0}/submit.py'.format(jobfolder)
        fp = open(fname,'w')
        fp.write('import os\n')
        fp.write('for i in range(0,{0}):\n'.format(counter))
        if self.JOBqueue=='SGE':
            cmd = '\'qsub jobscript{0}\'.format(i)'
            if simType=='transitions':
                cmd = '\'qsub -t 1-50:1 jobscript{0}\'.format(i)'
        elif self.JOBqueue=='SLURM':
            cmd = '\'sbatch jobscript{0}\'.format(i)'
        fp.write('    os.system({0})\n'.format(cmd))
        fp.close()

    def _extract_snapshots( self, eqpath, tipath ):
        tpr = '{0}/tpr.tpr'.format(eqpath)
        trr = '{0}/traj.trr'.format(eqpath)
        frame = '{0}/frame.gro'.format(tipath)
        
        gmx.trjconv(s=tpr,f=trr,o=frame, sep=True, ur='compact', pbc='mol', other_flags=' -b 2000')
        # move frame0.gro to frame50.gro
        cmd = 'mv {0}/frame0.gro {0}/frame100.gro'.format(tipath)
        os.system(cmd)
        
        self._clean_backup_files( tipath )
        
        
    def prepare_transitions( self, edges=None, bGenTpr=True ):
        print('---------------------')
        print('Preparing transitions')
        print('---------------------')
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            toppath = self._get_specific_path(edge=edge)
            
            for state in self.states:
                for r in range(1,self.replicas+1):
                    print('Preparing: PROT {0} {1} run{2}'.format(edge,state,r))
                    eqpath = self._get_specific_path(edge=edge,state=state,r=r,sim='eq')
                    tipath = simpath = self._get_specific_path(edge=edge,state=state,r=r,sim='transitions')
                    self._extract_snapshots( eqpath, tipath )
                    if bGenTpr==True:
                        for i in range(1,101):
                            self._prepare_single_tpr( tipath, toppath, state, simType='transitions',frameNum=i )
        print('DONE')  
        
        
    def _run_analysis_script( self, analysispath, stateApath, stateBpath, bVerbose=False ):
        fA = ' '.join( glob.glob('{0}/*xvg'.format(stateApath)) )
        fB = ' '.join( glob.glob('{0}/*xvg'.format(stateBpath)) )
        oA = '{0}/integ0.dat'.format(analysispath)
        oB = '{0}/integ1.dat'.format(analysispath)
        wplot = '{0}/wplot.png'.format(analysispath)
        o = '{0}/results.txt'.format(analysispath)

        cmd = 'pmx analyse -fA {0} -fB {1} -o {2} -oA {3} -oB {4} -w {5} -t {6} -b {7}'.format(\
                                                                            fA,fB,o,oA,oB,wplot,298,100) 
        os.system(cmd)

        
        if bVerbose==True:
            fp = open(o,'r')
            lines = fp.readlines()
            fp.close()
            bPrint = False
            for l in lines:
                if 'ANALYSIS' in l:
                    bPrint=True
                if bPrint==True:
                    print(l,end='')

    def run_analysis( self, edges=None, bProt=True, bParseOnly=False, bVerbose=False ):
        print('----------------')
        print('Running analysis')
        print('----------------')
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            print(edge)
            for r in range(1,self.replicas+1):

                # protein
                if bProt==True:

                    analysispath = '{0}/analyse{1}'.format(self._get_specific_path(edge=edge),r)
                    create_folder(analysispath)
                    stateApath = self._get_specific_path(edge=edge,state='stateA',r=r,sim='transitions')
                    stateBpath = self._get_specific_path(edge=edge,state='stateB',r=r,sim='transitions')
                    self._run_analysis_script( analysispath, stateApath, stateBpath, bVerbose=bVerbose )

        print('DONE')
