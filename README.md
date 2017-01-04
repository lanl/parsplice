The ParSplice code implements the Parallel Trajectory Splicing algorithm described in [1]. This method is part of the Accelerated Molecular Dynamics family of 
techniques developed in LANL over the last 16 years. These methods aim at generating high-quality trajectories of ensembles of atoms in materials. 
ParSplice uses multiple independent replicas of the system in order to parallelize the generation of such trajectories in the time domain, enabling 
simulations of systems of modest size over very long timescales. ParSplice includes capabilities to store configurations of the system, to generate and 
distribute tasks across a large number of processors, and to harvest the results of these tasks to generate long trajectories. ParSplice is a management 
layer that orchestrate large number of calculations, but it does not perform the actual molecular dynamics itself; this is done by external molecular dynamics 
engines.

This is a reference implementation of ParSplice that is not recommended for use in production. ParSplice is now under heavy development as part of the EXAALT 
project. Production quality EXAALT releases are expected to be available from Q4 2017. 


[1] Danny Perez, Ekin D Cubuk, Amos Waterland, Efthimios Kaxiras, Arthur F Voter, Long-time dynamics through parallel trajectory splicing, Journal of chemical theory and computation 12, 18 (2015)