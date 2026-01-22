# MatlabTgCalculation
TGCalc is a MATLAB script for simulating far-field patterns in transient grating spectroscopy.


### Description of files/directory structure ###

The source code consists of a collection of Matlab scripts (.m-files) and an accompanying input data files (provided as *.xlsx files).

The main directory contains the scripts of the core functionalities together with input data files that were used for initial silmulations. The files of the main directory are further described below. `TGCalc.m` is used to read Excel-based spreadsheet data and start calculation. 
`visualizeTGResult.m` provides output data visualization after calculation. Calculations are implemented in `TGCalc.m`. Note, all calculations may be done even without usage of `TGCalc.m` and `visualizeTGResult.m`. However, the usage of these files makes calculation process easier. 

Folder `inputs_example` contains examples of Excel spreadsheet used to check eligibility of the model.

Further scripts of core functionalities located in the main directory are the following:

`TGCalc.m`: main script used for reading input data files and forming MATLAB structures that are used as input parameters for simulations.

`processTGCalculation.m`: main program function, where program parameters are provided as inputs along with several supporting functions. Can be used separatelly by providing correct input data structures with command `some_result = processTGCalculation(grating, pump, pulse, sample, sys, io);`' in command window.
  `grating` structure contains following fields
    `slit`:       grating slit form factor and transmission type description
    `p`:          period of the grating parameter
    `df`:         grating duty factor
  `pump` and `probe` structures contain following fields
    `pd`:         grating phase depth
    `angle`:      incidence angle
    `wfc`:        wavefront curvature
    `waist`:      Gaussian beam inverted radius
    `lambda`:     wavelength
    `intensity`:  beam intensity
  `sample` structure contains following fields
    `t`:          sample thickness
    `n2`:         nonlinear complex refractive index
  `sys` structure contains following fields
    `g2s`:        grating to sample distance
    `ff`:         grating to PSD distance
    `psd`:        measurement aperture of PSD
    `ds`:         lateral coordinate separation
    `dsdp`:       division parameter for sys.ds
  `io` structure contains following fields
    `ddp`:        indicator used to save internal data
    `sdp`:        indicator used to save generated figures
    `wd`:         working path
    `filename`:   path to output file
    `i`:          additional information name parameter

visualizeTGResult.m: a script that includes sample routines for displaying the calculated results. Can be used separatly by running `visualizeTGResult(result);' command in command window after loading result of calculation *.mat file. Note, that result structure is created by TGCalc script.


### Installation instructions ### 

Included MATLAB scripts do not require any installtion, just copy to the working folder.


### Program execution ### 

To set up a simulation, an input TGCalc scripd must be started which allows to enter input data via GUI. The input parameters are read from the “Data” tab of the Excel file and are organized into six structlike columns: Grating (corresponding to the grating variable in the main function), Pump (pump), Probe (probe), Sample (sample), System (sys), and IO (io). At the start of program execution, the input data file filename.xlsx is
copied into the otput folder where calculation result will be saved as *.mat file also.


### Example of input data files ###

We provide users with a script that facilitates the definition of the required parameters to execute the main function, along with an example Excel data files containing a list of these parameters. These files contain the following examples: 

inputData\inputDataHeterodyne_dsdp.xlsx: Basic file used to estimate division parameter sys.dsdp

inputData\inputDataHeterodyne_tg-ff-g2s-<filename>.xlsx: Files used to simulate grating to sample distance effect on the far-field patterns

inputData\inputDataHeterodyne_tg-ff-pp<filename>.xlsx: Files used to simulate pump to probe ratio effect on far-field patterns formation

Full code description and usecases are provided within https://doi.org/10.1016/j.cpc.2025.109964.
