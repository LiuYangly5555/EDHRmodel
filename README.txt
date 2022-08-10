This repository contains the source code and data necessary to replicate the work in the manuscript " Evaporation-driven internal hydraulic redistribution alleviates root drought stress: mechanism and modeling".

---------------------
Parameter Description

In order to reduce the computational load of the model and the running time of the program (from about one hour to ten minutes), we scaled the plants and soil horizontally by a factor of 10: that is, the "plant rooting area" and "the vessel number" were divided by 10, and other parameters are the same as in the manuscript.

We have verified that horizontally scaling the plants and soil did not change the rate of soil drying, sap flow patterns, and water potential in each segment of the network, but only reduced the water flow rate in each segment of the model network by a factor of 10. Therefore, after exporting the data, all flow rates need to be multiplied by 10.

In addition, the units of the data exported by the program are not exactly the same as those in the manuscript and need to be converted. For example, the units of water potential and flow rate in progress are "cm" and "cm3 10min-1", they need to be converted to "kPa" and "cm3 h-1" as the manuscript.

---------------------
File Description

The files in the "Data" folder are the data used in the manuscript.

The files in the "Soil Function" and "Soil ET" folders are the functions for calculating soil characteristics and evaporation.

The files in the "Environmental Factors" folder are the set environment parameters including leaf water potential, relative humidity, wind speed, air temperature, soil water potential gradients.

The file in the "Network Flow" folder are the functions for solving the network flow rate in different conditions and hydraulic architectures.

---------------------
Program Running

The programs are divided into three categories. For convenience, the running of each program file is independent. All the above folders need to be loaded before running the program.

1. Simulate sap flow over time.
(In each program, soil parameters or plant parameters can be changed to simulate different scenarios）
Run "HA1.m" to simulate sap flow over time in HA1, the soil textures in shallow and deep layers are same.
Run "HA2.m" to simulate sap flow over time in HA2, the soil textures in shallow and deep layers are same.
Run "HA1_multi.m" to simulate sap flow over time in HA1, the soil textures in shallow and deep layers are different.
Run "HA2_multi.m" to simulate sap flow over time in HA2, the soil textures in shallow and deep layers are different.

2. Simulate instantaneous flow in different soil moisture.
Run "Instant_Flow_Soil.m" to simulate instantaneous flow rate in different shallow and deep soil moisture. By default, the program runs HA1. If you want to run HA2, you need to hide the code of HA1 and display the code of HA2. Besides, you can simulate instantaneous flow rate in different soil textures by changing soil parameters.

3. Simulate instantaneous flow in different plant factors.
Run "Instant_Flow_Hleaf.m" to simulate instantaneous flow rate in different leaf water potential.
Run "Instant_Flow_RLDs.m" to simulate instantaneous flow rate in different root length density of shallow layer.
Run "Instant_Flow_Lv.m" to simulate instantaneous flow rate in different vessel length.
Run "Instant_Flow_Droot.m" to simulate instantaneous flow rate in different vessel diameter of shallow root.
Run "Instant_Flow_rpit.m" to simulate instantaneous flow rate in different pit resistance.
Run "Instant_Flow_FcFp.m" to simulate instantaneous flow rate in different faction of the pit area.
