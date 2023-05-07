This project is the MATLAB code for all the experiments in the paper "6G Shared Base Station Planning Using an Evolutionary Bi level Multi objective Optimization Algorithm", consisting of three parts of experiments: 1) 6G Shared Base Station Simulation Experiment, 2) Model Sensitivity Analysis, and 3) Benchmark Testing Experiment.

1. 6G shared base station simulation experiment (folder TestonSimulation)
We have prepared a total of 5 test instances in the simulation experiment, SharedBS6G_ s1.m、SharedBS6G_ s2.m、SharedBS6G_ m1.m、SharedBS6G_ m2.m、SharedBS6G_ 2000.m, corresponding to the 5 test instances in Table 4 of the paper, each test instance was run 11 times.

Considering that calculating all test instances simultaneously takes a long time, in the code we provided, SharedBS6G_s1.m is calculated by default, and this test instance is only run once by default. You can determine the number of runs of the test instance by setting the value of the parameter Run_Times in the TESTonSimulation\ExperimentalDesign_SBS.m file, such as Run_ Times = 11。

If you would like to verify other test instances, you can determine which test instances to run by setting the value of the parameter Problems in the TESTonSimulation\ExperimentalDesign_SBS.m file.
For example:
Problems = {...   
'SharedBS6G_s1'...
'SharedBS6G_s2',...
'SharedBS6G_m1',...
'SharedBS6G_m2',...
'SharedBS6G_2000',...
}

We also provide a simple test instance SharedBS6G.m for program debugging. The Problems array in the TESTonSimulation\ExperimentalDesign_SBS.m file only includes this test instance, which can quickly debug the program.
For example:
Problems = {...   
'SharedBS6G'
}

Note that after setting the test instance, add the TestonSimulation folder and all its subfolders to the MATLAB path, and then run the TESTonSimulation\ExperimentalDesign_SBS.m file to start this experiment.

After the program runs, the detailed data will be saved in the TESTonSimulation\PlatEMO\Data\stMOBEA folder, with the following data structure:

result

|-Number of the upper level evaluations

|-Number of the lower level evaluations

|-

​   |--Number of the upper and lower level evaluations for each round of evolution

​   |--Solutions obtained by algorithm

​   |--IGD

​   |--HV

​   |--Non-dominated solutions for each round of evolution

|-run time


2. Model Sensitivity Analysis (folder TestonSimulation)

We have prepared a total of 9 test instances SharedBS6G_t1.m、SharedBS6G_t2.m、SharedBS6G_t3.m、SharedBS6G_t4.m、SharedBS6G_t5.m、SharedBS6G_t6.m、SharedBS6G_t7.m、SharedBS6G_t8.m、SharedBS6G_t9.m, divide each of the three test instances into a set of experiments, corresponding to the three sets of experiments in Table 2 of the paper. If you would like to verify these instances, please set the relevant parameters in the TESTonSimulation\run.m file.

For example (conducting the first set of experiments):

main('-algorithm',@TEST_123,'-problem',@SharedBS6G_t1,'-Nu',1,'-Nl',1,'-run',1,'-save',1);

main('-algorithm',@TEST_123,'-problem',@SharedBS6G_t2,'-Nu',1,'-Nl',1,'-run',1,'-save',1);

main('-algorithm',@TEST_123,'-problem',@SharedBS6G_t3,'-Nu',1,'-Nl',1,'-run',1,'-save',1);

For example (conducting a second set of experiments):
main('-algorithm',@TEST_456,'-problem',@SharedBS6G_t4,'-Nu',1,'-Nl',1,'-run',1,'-save',1);

main('-algorithm',@TEST_456,'-problem',@SharedBS6G_t5,'-Nu',1,'-Nl',1,'-run',1,'-save',1);

main('-algorithm',@TEST_456,'-problem',@SharedBS6G_t6,'-Nu',1,'-Nl',1,'-run',1,'-save',1);

After setting the parameters, please run the TESTonSimulation\run.m. After the program runs, the detailed data will be saved in the TESTonSimulation\PlatEMO\Data\stMOBEA folder.

3. Benchmark test experiment (folder TestonBenchmark)

The benchmark testing experiment is to test the performance of the proposed SABLEA-PM and other comparative algorithms on benchmark problems (TP1-2 and DS1-5), corresponding to the experiment in section 5.3 of the paper. Each algorithm was run 21 times on each benchmark problem.

Run the TESTonBenchmark\ExperimentalDesign.m file to execute the benchmark test experiment. After the program runs, the detailed data of each algorithm will be saved in the corresponding folder located in TESTonBenchmark\PlatEMO-bilevel\Data, and its data structure is the same as that in the 6G shared base station simulation experiment.