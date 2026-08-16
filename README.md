# Friction biofilms

A C++ simulation framework for modeling the impact of friction in bacterial colonies/biofilms.
Project done by Lucca Johann Leal for a MSc in Theoretical Physics dissertation for the University of Edinburgh (2025-2026)

This code is built upon the framework designed by Rory Claydon (https://github.com/roryclaydon1994/BiofilmDES?tab=readme-ov-file#references) and Laila Saliekh (https://github.com/lailasaliekh/Candida-project).




## Table of Contents

- [About](#about)
- [Building](#building)  
- [Running Simulations](#running-simulations)    
- [Slurm Scripts](#slurm-scripts)
- [Parameters](#parameters)
- [Analysis](#analysis)
---

## About

A Discrete Element Simulation to study the impact of friction in growing bacterial coloniescolonies
---

## Building
To compile the code run from the terminal (Linux operating system) one needs CMake (v >= 3.16.0) for compilation and assumes g++ 9.3.0 or above.
Go to the main directory "./Friction_Biofilms" and use the command
```bash

./quickMake.sh
```
Remember to recompile once changing the hard-coded parameters

## Running Simulation
From the terminal inside the main directory, use:
```bash
./build/Main/main.out <output_dir>
```

## Slurm Scripts
In order to submit jobs to the slurm cluster we can either use the following, 
for a single parameter set submission
``` bash
sbatch one_job.sh
```
For loop submitting multiple repeats
``` bash
./submit_repeats.sh <number_of_repeats>
```
To allow these to be run as a program or script we use,
``` chmod +x <script.sh> ```

## Parameters
The initial number of cellsand the parameters of each type can be changed from the ```main.cpp``` file
For the number of cells
``` C++
  int numTypeA = 1;      
  int numTypeB = 1;     
```


 ```main.cpp```
``` C++
            if (isTypeA) {
                auto* rod = new RodShapedBacterium{
                    x, y, 0,  //random position x, y, z (in 2D z=0)
                    angle,// Random angle
                    constants::pi * 0.5, 
                    RodShapedBacterium::mAvgGrwthRate, // Type A growth rate
                    4, // Type A initial length
                    0, // non-chaining, 1 for chaining
                    0.5, //radius of the cell
                    constants::Lambda1 // Lambda 1 for standard drag
                };
                initial_conditions.push_back(rod);
                numTypeA--;
            } 
            else if (isTypeB) {
                auto* rod = new RodShapedBacterium{
                    x, y, 0, //random position x, y, z (in 2D z=0)
                    angle,// Random angle
                    constants::pi * 0.5, 
                    RodShapedBacterium::mAvgGrwthRate, // Type B growth rate
                    4, // Type B initial length
                    0, // non-chaining, 1 for chaining
                    0.5, //radius of the cell
                    constants::Lambda2 // Lambda > 1 for higher drag , for base drag = 1
                };
                initial_conditions.push_back(rod);
                numTypeB--;
            }
```

To change the parameters of the bacteria, such as the growth rate, open the ```constants.hpp``` file inside the includes directory.

## Analysis
The analysis aspect of the project is contained within the Analysis directory. In this directory, any files not mentioned in this read me have been previously developed by the previous builders of the code (see first line).

The data analysis is split into several files that inherit functions from one another. The most basic one is ```drag_functions.py``` which has standard calculations and methods used. Then, ```plotting_functions.py``` inherits various functions and plots the different variables in a matplotlib figure. ```Automatic_plotting.py``` and ```Plot_channels.py``` call the functions from ```plotting_functions.py``` and organise the data whilst setting labels.

Then, to visualise the data, jupyter notebooks (```Channel_analysis.ipynb```, ```FreeGrow_plots.ipynb``` and ```Interaction_plots.ipynb```) call the different functions to produce the data visualisation. 

There are some other files that have been used to create graphs specifically for the analytical model of my disseration, those are ```Analytical_wrapping.ipynb```, ```analytical_model_functions.py``` and ```report_figures.py```. These follow a similar structure but need some reorganisation.