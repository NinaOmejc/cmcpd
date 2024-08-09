-----

First run preprocessing files inside "preprocessing" folder. 
The scripts should be ran as they are in the alphabetic order. 

-----

Then, for connectivity:

1. Open "create_config.m" script, adjust root_path and run it (change other settings if required, e.g. which tasks to run)

2. Config file will be created and saved inside ".\configs" folder.

3. Run "main_connectivity.m" script with two arguments: "config_path" and "job_idx", e.g.
    
	main_connectivity('D:\Experiments\corticomuscular_analysis\configs\config_loc_test_cwt_pc_dbi.mat', 1)
	
	The results are saved inside folder ".\data\real\" 

	If you want to run a method for both subjects and all tasks, either use job_idx = 0 (it's then run in parallel using matlab parallel
	toolbox, if you have) or run jobs in for loop:
	
	for ijob = 1:8  % ten because 2 subjects and 4 tasks for each (no C task included)
		main_connectivity(<config_path>, ijob)
	end
	

4. For analysis of results see and run individual scripts inside '.\code\analysis' folder. 
    "<method>_group_analysis" scripts collect results from the individuals and join them in a single struct.
    "<method>_group_plots" scripts load single group struct and make plots.
    Note, the path to config file inside the scripts should be changed, because the settings are imported as well.
    The results are saved inside folder ".\analysis\<method>\"
	
	Currently, there are the data of two subjects, one healthy and one patient, because all "analysis" scripts can only be 
	ran if you have both groups. 
	
	I commented out all plotting functions that dont work for the particular configs (topo plots are not possible with only 
	one channel, avgcen, and boxplots are not possible with only one value per group)
	
	Also, to run the group plotting scripts, the EEGLAB toolbox is needed (https://sccn.ucsd.edu/eeglab/download.php)
	
	----
	
	for more info, feel free to ask me via email: nina.omejc@ijs.si