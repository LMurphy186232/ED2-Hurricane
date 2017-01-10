#!/bin/bash

#----- Main settings.  Do look at run_sitter.sh and epost.sh for additional settings. -----#
here=$(pwd)            # Current path.
email="none"           # Your e-mail. Put a backslash before the @
queue="anyqueue"       # Queue to run run_sitter.sh and epost.sh
runtime="infinite"     # Run time request
memory=2048            # Requested memory (Mb)
sbatch=$(which sbatch) # SLURM command to submit job.
rscript="nothing"      # Which script to run with epost.sh.  See options below. 
                       #     Multiple scripts are allowed, put spaces between them.
                       #     (e.g. rscript="plot_monthly.r plot_fast.r")
#------------------------------------------------------------------------------------------#
#     Which scripts to run.                                                                #
#                                                                                          #
#   - read_monthly.r - This reads the monthly mean files (results can then be used for     #
#                      plot_monthly.r, plot_yearly.r, and others, but it doesn't plot any- #
#                      thing.)                                                             #
#   - yearly_ascii.r - This creates three ascii (csv) files with annual averages of        #
#                      various variables.  It doesn't have all possible variables as it is #
#                      intended to simplify the output for learning purposes.              #
#   - plot_monthly.r - This creates several plots based on the monthly mean output.        #
#   - plot_yearly.r  - This creates plots with year time series.                           #
#   - plot_ycomp.r   - This creates yearly comparisons based on the monthly mean output.   #
#   - plot_povray.r  - This creates yearly plots of the polygon using POV-Ray.             #
#   - plot_rk4.r     - This creates plots from the detailed output for Runge-Kutta.        #
#                      (patch-level only).                                                 #
#   - plot_photo.r   - This creates plots from the detailed output for Farquhar-Leuning.   #
#   - plot_rk4pc.r   - This creates plots from the detailed output for Runge-Kutta.        #
#                      (patch- and cohort-level).                                          #
#   - plot_budget.r  - This creates plots from the detailed budget for Runge-Kutta.        #
#                      (patch-level only).                                                 #
#   - plot_eval_ed.r - This creates plots comparing model with eddy flux observations.     #
#   - plot_census.r  - This creates plots comparing model with biometric data.             #
#   - whichrun.r     - This checks the run status.                                         #
#                                                                                          #
#   The following scripts should work too, but I haven't tested them.                      #
#   - plot_daily.r   - This creates plots from the daily mean output.                      #
#   - plot_fast.r    - This creates plots from the analysis files.                         #
#   - patchprops.r   - This creates simple plots showing the patch structure.              #
#   - reject_ed.r    - This tracks the number of steps that were rejected, and what caused #
#                      the step to be rejected.                                            #
#------------------------------------------------------------------------------------------#


#----- Make sure e-mail and queue are pre-defined. ----------------------------------------#
if [ "x${email}" == "xnone" ]
then
   echo " You must set up variable \"email\" before running sim_sitter.sh"
   exit 99
fi
if [ "x${queue}" == "xanyqueue" ]
then
   echo " You must set up variable \"queue\" before running sim_sitter.sh"
   exit 99
fi
if [ "x${rscript}" == "xnothing" ]
then
   echo " You must set up variable \"rscript\" before running sim_sitter.sh"
   exit 99
fi
#------------------------------------------------------------------------------------------#



#----- Make substitutions. ----------------------------------------------------------------#
sed -i~ s@"here=\"/nowhere\""@"here=\"${here}\""@g          ${here}/run_sitter.sh
sed -i~ s@"recipient=\"none\""@"recipient=\"${email}\""@g   ${here}/run_sitter.sh
sed -i~ s@"here=\"/nowhere\""@"here=\"${here}\""@g          ${here}/epost.sh
sed -i~ s@"queue=\"anyqueue\""@"queue=\"${queue}\""@g       ${here}/epost.sh
sed -i~ s@"rscript=\"nothing\""@"rscript=\"${rscript}\""@g  ${here}/epost.sh
#------------------------------------------------------------------------------------------#


#----- Job preferences. -------------------------------------------------------------------#
joblog="${here}/out_sitter.out"
jobpref=$(basename ${here})
jobname="${jobpref}-sitter"
#------------------------------------------------------------------------------------------#



#----- Submit run_sitter.sh in batch mode. ------------------------------------------------#
${sbatch} -p ${queue} --mem-per-cpu=${memory} -t ${runtime} -o ${joblog} -J ${jobname}     \
   -n 1 --wrap="${here}/run_sitter.sh"
#------------------------------------------------------------------------------------------#
