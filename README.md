# MKJob_matlab
Generate job files for a Matlab program on Argon cluster


## v1.0 (10/31/2019)

Here is an example about how to generate, submit and manage a serier of Matlab jobs.

* **run_angles.m** and **gd_function.m**are the simulation files.

* **create_matlab_jobs.py** is saved in the same folder with the m files. This command generates many subfolders to cover a series of simulaitons.

* **dwt_matlab-job-file.job** is the template for the job files.

* **modify-matlab-jobfile** is a command which can modify the resource request in the job file.
