This code is free to use for research purpose, and available via the GNU Public License. If you find this useful in your research, please cite the publication listed below.

Wasseem Al-Rousan, Caisheng Wang, and Feng Lin,"Modular Control of Discrete Event System for Modeling and Mitigating Power System Cascading Failures". arXiv, 2025. [https://doi.org/10.48550/arXiv.2504.07496]
 

# Modular_DES
Mitigation of power system cascading failure using modular discrtete event systems
There are two main parts for this work: 

#Part 1: 
Preperation of supervisors based on the modular control approach of the power system using aut_files folder
    the files required to prepare the supervisors can be found under the aut_files folder
    First the power system components are translated into automatas , then parallel operation , specification then supervisors synthesis operations are performed
    these operations are performed by using LibFAUDES C++ library (more info can be found on  [https://fgdes.tf.fau.de/faudes/faudes_about.html]  and   
    [https://ieeexplore.ieee.org/document/4605933] )
    These operations are performed using the matlab files that use LibFAUDES library using MEX, an example of using mex command
    
                      mex -lfaudes -L"D:\IEEE_300_Bus\libFAUDES" Automaton.cpp
                      
    which is located in createautomatas.m file, after Automaton.cpp is compiled into mex format , it can be used iteratively as a function to generate Automatas for the all      the power system components as illstrated in the paper for all the lines, generators and loads in the given case study of the power system. The Automatas of all buses in     the systems and associated components with each bus are preapred using he file : prepare_plants_w_paralle_newl.m. In similar way, specifications Auotomatas and     Supervisors Automatas are prepared. 
    
    
    After supervisors are synthesized, another Matlab file (Supervisors_Matlab_info_general.m) transform the supervisors to directed graphs format in Matlab. 
    Note: faudes.dll file need to be copied to the aut_files folders to enable the use of libfaudes library. also you may need to modify the header files directory in the        cpp files. 
The "off-line" operations on DES can be summrized in Fig. 1 below

<img width="150" height="250" alt="flowchart_1" src="https://github.com/user-attachments/assets/e5c874d8-515f-43b3-a4b8-e9cf51a4cd2b" />

Fig. 1. Flowchart of the supervisor synthesis procedure

#Part 2: 
employing the synthsized superviros into the smulations of power systems cascading failure using sim_files folder
  the file DESsimuatoroverall_DES_Modular_MC_100.m peform Monte Carlo (MC) simulations of the casacding failure process with and without control 
  one of the control methods used is the DES modular control:
  After preparing teh supervisors in Part 1 above, the superviros are used using supervisors iterators files in the modular control approach. 
  The casacding failure simualtion is absed on the DCSIMSEP and ACSIMSEP algorithms can be found on [https://github.com/phines/acsimsep] 
  changes were added in the distibuted control file where we added our approach

  The overall co-simulation of both DES and power systems can be summrized in Fig. 2 below. 

  <img width="200" height="400" alt="flowchart_2" src="https://github.com/user-attachments/assets/a6b8576e-1144-4299-bdb3-65f2ecadcc3e" />
  
  Fig. 2. Flowchart of the implementation of the modular DES controller

