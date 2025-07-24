I will develop the read me more later this week but for now, here are the steps to run oscar. 

1. Organize your datasets and diagnostics folders. The diagnostics folder has the checks along the way to see the interpolated graphs and all the gradients, I will work on creating a flag to turn this off. These folders should be on the same level as pyoscar. For example
    |_datasets
    |_diagnostics
    |_pyoscar

    a. How to organize the datasets folder
        |_datasets
            |_CURRENTS
                |_FINAL (EMPTY RN)
                |_INTERIM
                    |_PODAAC
                        |_CMEMS
                            |_2020
                                |_01 
                                    |_oscar_currents_nrt{year}{month}{day}
                        |_NEUROST
                            |_2020
                                |_01
                                    |_oscar_currents_nrt{year}{month}{day}
            |_SSH
                |_CMEMSNRT
                    |_SRC
                        |_2020
                            |_01
                                |_nrt_global_allsat_phy_l4_{year}{month}{day}_*.nc
                |_NEUROST
                    |_SRC
                        |_2020
                            |_01
                                |_NeurOST_SSH-SST_{year}{month}{day}_*.nc
            |_SST
                |_SRC
                    |_CMC
                        |_2020
                            |_01
                                |_{year}{month}120000-CMC-L4_GHRSST-SSTfnd-CMC0.1deg-GLOB-v02.0-fv03.0.nc
            |_WIND
                |_ERA5
                    |_SRC
                        |_2020
                            |_ERA5SRC202001.nc

        b. How to organie your diagnostics folder

            |_diagnostics
                |_currents
                    |_figures
                        |_2020
                            |_earth (will change to globe lol)
                                |_cmems
                                    |_
                                |_neurost
                                    |_
                            |_gulfstream
                                |_...
                            |_...
                |_final
                |_interim
                |_nrt 
                    |_interp_grad
                        |_figures
                            |_the sub folders should populate automatically


2. pyoscar/config/io_config and fill in your paths

3. cd tp the folder outside of pyoscar and run this in your terminal >> python3 -m pyoscar.main
    To note: if you want to plot the currents run >> python3 -m pyoscar.analysis.main
    (I havent made this automatic yet and the analyis folder is not finished developing as yet so I have just been observing the final plots this way. Also, my folder oaths are hardcoded in here for right now, you will see in analysis/main.py - if you want to plot your currents, you can remove my hardcoded paths for right now )


More about how pyoscar is organized right now

pyoscar
    |_analysis (this is where I plan to do all the analyses including plotting the currents, validation agains the drifter etc ... I have the drifter code for the most part but gotta integrate it)
    |_computation
        |_ compute_currents.py (this is like the "main" for the computation, a short function encapsulating everything)
        |_constants.py (all the constants)
        |_equatorial_treatment.py(all the functions for the equatorial treatment)
        |_forcing.py(all the functions for the forcings)
        |_interp_and_grad.py(where we interpolate the the ds and calculate the gradient)
        |_physics.py(all the other functions)
    |_config
        |_io_config.yaml (this is the only file that needs to be adjusted for the currents)
        |_setup.py(where I turned the yaml stuff into python variables)
    |_test (I have some preliminary stuff here but nothing major yet)
    |_utils
        |_data_utils.py (where all the data/file operations are)
    |_main.py (the entry point)
    |_README.md
 ps: you will notice that all the folders and sub olders have __init__.py,these files are empty and serve the purpose of marking each directory as a Python package, allowing for proper module imports throughout the project.


