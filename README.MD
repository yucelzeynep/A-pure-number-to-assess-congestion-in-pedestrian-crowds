There are two versions of the cn_computation code, cn_computation and cn_computation2.

cn_computation computes cn, rotor, density etc. only in a given area chosen by the user through the file parameters, while cn_computation2 computes it everywhere (in a maximum rectangular area containing all data) but also performs some statistics limited to an area chosen by the user (through parameters2).

The data are in the directory positions, and correspond to the "low density realistic setting" in the PRX manuscript. There are 10 repetitions of the same experiment, and the parameter ripet=10 is passed through the parameters files.

DPIX and CELL_PIX give the size in pixels of a cell in the figures.

GRID_SIZE is the actual size of a grid cell, TIME_INTERVAL is the time step over which the integration of the velocity field is performed, CN_COMP_RAD decides the size of the ROI (it's the radius in GRID_SIZE units). X_MIN, Y_MIN etc. in parameters give the area in which data is checked, and TI the time span (end of experiment, beginning is supposed to be 0), These are not needed in parameters2 (that can analyse also negative time), but we can pass an area of interest with XI_MIN etc.
Analysed data is written in the folders data and data2
