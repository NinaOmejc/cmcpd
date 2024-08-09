Preprocessing part 

First run preprocessing files inside "preprocessing" folder. 
The scripts should be ran as they are in the alphabetic order. 


Connectivity part

1. Open "create_config.m" script, adjust root_path and run it (change other settings if required, e.g. which tasks to run)

2. Config file will be created and saved inside ".\configs" folder.

3. Run "main_connectivity.m" script with two arguments: "config_path" and "job_idx", e.g. 
main_connectivity('D:\Experiments\corticomuscular_analysis\configs\config_loc_test_cwt_pc_dbi.mat', 1)	
The results are saved inside folder ".\data\".

For more info, feel free to ask me via email: nina.omejc@ijs.si
