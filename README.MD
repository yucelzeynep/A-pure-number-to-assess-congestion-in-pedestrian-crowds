
# A pure number to assess “congestion” in pedestrian crowds

This repository contains the resources necessary for reproducing the results reported in our manuscript submitted to PRX (2020). 

A preprint of the -earlier- version of the article can be found at [Zanlungo2020], and the interested reader may also refer to [Feliciani2018] for the sprouting of the idea.

The implementation in this repository is the original work of Dr. Francesco Zanlungo. Any comments or questions can be communicated to zanlungo AT_MARK atr.jp.


**Main routines**

The programs are developed in C++ with no specific dependencies. But Simple DirectMedia Layer (SDL1.2) is necessary for producing the graphics. 

There are two versions of the cn_computation code, cn_computation and cn_computation2.

cn_computation computes cn, rotor, density, etc. only in a given area chosen by the user through the file parameters, while cn_computation2 computes it everywhere (in a maximum rectangular area containing all data) but also performs some statistics limited to an area chosen by the user (defined in the file parameters2).

**Hyperparameters**

Hyperparameters are defined in the files parameters and parameters2 (for cn_computation and cn_computation2, respectively)

DPIX and CELL_PIX give the size in pixels of a cell in the figures.

GRID_SIZE is the actual size of a grid cell.

TIME_INTERVAL is the time step, over which the integration of the velocity field is performed.

CN_COMP_RAD is is the radius in GRID_SIZE units and it decides the size of the ROI (Region Of Interest).

X_MIN, Y_MIN etc. in parameters give the area, in which data is checked.

TI is the time span (end of experiment, the beginning is supposed to be 0), These are not needed in parameters2 (that can analyse also negative time), but we can pass an area of interest with XI_MIN etc. 

**Output data**

Analysed data is written into the folders /data and /data2.
The data in the folder /positions correspond to the "low density realistic setting" in the manuscript. There are 10 repetitions of the same experiment, and the parameter ripet=10 is passed through the parameters files (to let the program know how many repetitions are in the data set).


**References**

[Zanlungo2020] Zanlungo, F., Feliciani, C., Yucel, Z., Nishinari, K., & Kanda, T. (2020). Some considerations on crowd Congestion Level. arXiv preprint arXiv:2004.01883.

[Feliciani2018] Feliciani, C., & Nishinari, K. (2018). Measurement of congestion and intrinsic risk in pedestrian crowds. Transportation Research Part C: Emerging Technologies, 91, 124-155.
