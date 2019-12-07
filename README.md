# ROBUST DIGITAL FILTERS MINIMAX FIR FILTER

Implemented following paper:
Robust Digital Filters Part 1 — Minimax FIR Filters

https://www.researchgate.net/publication/332815480_Robust_Digital_Filters_Part_1_-_Minimax_FIR_Filters

This project implements a robust digital filter, where robust digital filters are the filters that offers optimal performance under variation of filters parameters. 
The characteristics performance of robust digital filters against the parameter uncertainties are discussed in project report. 
The performance and design formulations of robust FIR filters in L8 (minimax) are discussed as nonsmoothed convex problems. 
In this project, L8 -robust FIR filter is designed (optimized) using an accelerated sub-gradient algorithm and implemented in MATLAB  


More details in project report.pdf

## Getting Started

The code for the project implemented is in file "Final.mat".
 

Enter the values for the following varible before running the main program, default values are given in program:

%fp:	passband edge 

%fa:	stopband edge

%N:	 filter length 

%gam:	gama variable

%b:	Beta variable

%M:	Number of frequency grids

%K:	number of iteration for Size of bounding box


## References

[1] Wu-Sheng Lu ; Takao Hinamoto , "https://www.researchgate.net/publication/332815480_Robust_Digital_Filters_Part_1_-_Minimax_FIR_Filters".
doi: 10.1109/LSP.2013.2267593

[2] W.-S. Lu, Course notes of Advanced Mathematical optimizations. 




