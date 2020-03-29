How to use micro_model developed by Nid

****************************** Focus on micro_model_3D only ******************************
- micro_3D_uniaxial & micro_3D_submodel = master script
- model_example_3D = all related files to show process
___________________________________________________________________________________________

micro_3D_uniaxial

1. 01_modelGeom, 02_modelSetting_Mises, DP1000M, surface_file: to generate .cae and .inp
2. Save .cae for further use
3. Submit simulation .inp for output .odb
4. 03_postProc_all: on .odb to get variables in excel lode, mises, peeq, triax, volume
5. getLocalElems: on .cae to get excel localElems (only one point canbe selected)
6. HET_analysis_uniaxial.ipynb: process all data for final result
>> input files <<
- triax, lode, peeq, volume: variables
- localElems: labels of local elements
- compPEEQ_outPunch: strain path of macro model for comparison
>> output files <<
- global-localVar: averaged variables of global and local elements
- PEEQ-displace: progress of PEEQ against boxsize displacement
- surfFac: surface factor result

___________________________________________________________________________________________

micro_3D_submodel

1. component_model.odb: as reference to submodel
2. 01_modelGeom, 02_modelSetting_Mises, DP1000M, surface_file: to generate .cae and .inp
3. Save .cae for further use
4. Submit simulation .inp for output .odb
5. 03_postProc_all: on .odb to get variables in excel lode, mises, peeq, triax, volume
6. getLocalElems_pntX: on .cae to get excel localElems_pntX (multiple local points canbe selected)
7. HET_analysis_submodel.ipynb: process all data for final result
>> input files <<
- triax, lode, peeq, volume: variables
- localElems_pntX: labels of local elements
- compPEEQ_Xum_middleStrip: strain path of macro model for comparison (30um for roughness, 450um for waviness)
>> output files <<
- globalVar: averaged variables of global elements
- localVar_pntX: averaged variables of local elements
- surfFac_allpnts: globalPEEQ / localPEEQ at every timesteps
- surfFac: final result of surface factor
- DIL-3D: variables and DIL plot in 3D

___________________________________________________________________________________________