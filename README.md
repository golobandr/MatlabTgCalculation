# MatlabTgCalculation
Computation of diffraction patterns of light in a transient grating (TG) geometry scheme. The output intensity distribution is calculated based on the diffraction integral in the Fresnel and Fraunhofer regimes.

The code is designed to assist both experimental planning and post-processing interpretation by modeling the optical response of TG configurations across a wide range of conditions. It supports input through structured MATLAB variables or Excel-based spreadsheets and provides automated consistency checks and visual output generation. The implementation includes integration over detector pixels, enabling realistic simulations that account for spatial averaging and resolution effects. We demonstrate the software’s capabilities through representative use cases, including the influence of the grating-to-sample distance, the pump-to-probe intensity ratio, and the selection of the division parameter governing pixel integration accuracy. The code is freely available and modular, facilitating its adaptation to different experimental geometries and beam conditions. While full validation is provided elsewhere, this work establishes the core methodology and illustrates the practical value of the tool for TG spectroscopy research.

`TGCalc.m` is used to read Excel-based spreadsheet data and start calculation. 
`visualizeTGResult.m` provides output data visualization after calculation. 
Calculations are implemented in `TGCalc.m`. Note, all calculations may be done even without usage of `TGCalc.m` and `visualizeTGResult.m`. However, the usage of these files makes calculation process easier. 

Folder `inputs_example` contains examples of Excel spreadsheet used to check eligibility of the model.


Full code description and usecases are provided within https://doi.org/10.1016/j.cpc.2025.109964.
