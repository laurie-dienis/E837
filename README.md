# Analysis of experiment e837
Code of the analysis of the experiment e837, a resonant elastic scattering experiment using ACTAR TPC. The code is using [ACTROOT](https://github.com/loopset/ActRootV2.git) to run.

## To convert the data : 

```
actroot -r tpc
```
```
actroot -r sil
```
```
actroot -f
```
```
actroot -m
```

## To visualize the data event by event: 
```
actplot -v
```

## To lauch the post-analysis code : 
```
cd PostAnalysis/
root 'Runner.cxx("number_of_the_pipeline")'
```
Pipepline 1 : PID particle identification, DeltaE-E  
Pipeline 2 : Computing the excitation energy of the compound nucleus   
Pipeline 3 : Plotting the Kinematics energy vs angles  
Pipeline 4 : NeuralNetwork form carbon induced scattering 
