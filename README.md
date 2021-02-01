# Ksp_phase_separation
<h2> 1. NFsim simulations </h2>
<ul> 
  <li> We generate the <b> two components </b> bngl files (http://bionetgen.org/) using a custom python script (file 00); <b>three components</b> files (A-B-C system) are written manually.  </li>  
 <li> The bngl file (file 01) is then converted to an xml file (file 02); detailed methods can be found here https://emonet.biology.yale.edu/sites/default/files/files/NFsim_manual_v1_11.pdf.  </li>
   <li> The xml file is then used as an input file to run multiple stochastic trajectories. The jobs are parallely executed in a high performance computing facility (https://health.uconn.edu/high-performance-computing/) using a custom shell script (file 03). </li>
    <li> Each trajectory or run would generate two output files: <ul> <li> gdat file containing the observables data (Run_x.gdat) </li> <li> species file containing the molecular clusters at the last time point of simulations (Run_x.species)</li> </ul> </li>
  <li> Once we have multiple stochastic runs, we use another python script (file 04) to compute two variables: <ol> <li> average counts of obvervables (Free or Total molecular concentration) across multiple runs </li> <li> distributions of molecular clusters at the last time point (note that the distribution does not correspond to a steady state when total concentration diverges with time in "CMC" simulations) </li> </ol> </li> 
  <li> For <b>"FTC" method</b>, we gradually increase molecular counts (100, 200, 300, ..., 1000) at fixed Kd (3500 molecules) to titrate up the concentrations; creation and decay rates are set to zero to make it a closed system. </li>
  <li> On the other hand, for <b> "CMC" method</b>, we begin with 1 molecules each (Kd = 3500 molecules) and increase the creation rates (1290000, 1300000, 1310000 molecules/s) with a fixed decay rate (10000 1/s) to titrate up the clamped concentrations. File_05 is an example for CMC model for 4v-4v system. </li>
  <li> <b>MixedValent system:</b> File_06 and File_07 are examples for 5v-3v and A3 - B13 - C6 systems respectively.  </li>
</ul>
<h2> 2. SpringSaLaD simulations </h2>
<ul>
    <li> SpringSaLaD models are built using the graphical user interface (GUI). The jar file and detailed tutorail can be found here https://vcell.org/ssalad. </li>
    <li> After defining the molecules and reaction rules, volume of the simulation box and simulation time (along with timestep and output inteval time) are set. File_08 is an example input file (for the "reference system") containing all the simulation details. </li>
  <li> SpringSaLaD has funtionalities to launch the multiple simulations locally from the GUI. But due to the intensive computations, we run multiple simulations (file_09) in our hpc facility. File_10 is the numerical solver needed to execute such simulations. </li>
  <li> Once multiple runs are executed, a python script is used to compile all those data and perform statistical analysis of the molecular concentrations (file_11) and molecular clusters (file_12).</li>
    
      
