# Description
### This is a script generating the input file "wt.in" of WannierTools.
### Put this script into the folder of constucting the Wannier functions and command:
python gen_wt_in.py
### Notice, before using, you should set the parameters in the beginning of the script, especially the "calc_task" and "wt_path".
### If you want to set "calc_task=24", i.e. calculate the CPGE tensor, because the origin WannierTools never have this function, I recommand you install my script "CPGE.f90" into WannierTools:
