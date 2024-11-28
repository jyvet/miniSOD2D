Splitting full_convec_ijk into 3 subroutines in order to reduce register spilling.


     Number of elements          =  100000
     Element order               =  4
     Number of nodes per element =  125
     Nodes on mesh               =  12500000
     Number of runs              =  100


     Timings:
     ----------------------------------------
     Avg. convective time     =  9.76052755499790675E-2
     Avg. TET convective time =  0.52057270731000704
     Avg. diffusive time      =  2.14014072899885825E-2
     Max. convective time     =  0.89435866300004818
     Max. TET convective time =  0.52119506099984392
     Max. diffusive time      =  2.14360630000101082E-2
     Min. convective time     =  8.78463849999207014E-2
     Min. TET convective time =  0.47557626299999356
     Min. diffusive time      =  2.13360509999347414E-2
     ----------------------------------------
     Variation convec.        =  8.2629988333689042
     Variation TET convec.    =  8.76319433563462619E-2
     Variation diffu.         =  4.67315063538609795E-3
     ----------------------------------------

     Basic results:
     ----------------------------------------
     --| Convec
     Max Rmass     =  4.057812691E-2 Min Rmass     =  -1.48370039
     Max Rmom(:,1) =  0.203694701 Min Rmom(:,1) =  -9.26724052
     Max Rmom(:,2) =  0.203694716 Min Rmom(:,2) =  -9.26724052
     Max Rmom(:,3) =  0.203694671 Min Rmom(:,3) =  -9.26724052
     Max Rener     =  14.5912218 Min Rener     =  -313.078461
     ----------------------------------------
     --| Diffu
     Max Dmass     =  0.122479714 Min Dmass     =  -0.464846373
     Max Dmom(:,1) =  17.5130577 Min Dmom(:,1) =  -10.0694799
     Max Dmom(:,2) =  17.5131435 Min Dmom(:,2) =  -10.069479
     Max Dmom(:,3) =  17.5126896 Min Dmom(:,3) =  -10.0697556
     Max Dener     =  518.689026 Min Dener     =  -225.701828
     Basic TET results:
     ----------------------------------------
     --| Convec
     Max Rmass     =  27.4684582 Min Rmass     =  -24.6873055
     Max Rmom(:,1) =  210.747665 Min Rmom(:,1) =  -195.865738
     Max Rmom(:,2) =  204.161804 Min Rmom(:,2) =  -196.02121
     Max Rmom(:,3) =  204.975861 Min Rmom(:,3) =  -197.575928
     Max Rener     =  67.3746109 Min Rener     =  -116.874619
     ----------------------------------------
     --| Diffu
     Max Dmass     =  -0.372548997 Min Dmass     =  -0.372548997
     Max Dmom(:,1) =  -0.372548997 Min Dmom(:,1) =  -0.372548997
     Max Dmom(:,2) =  -0.372548997 Min Dmom(:,2) =  -0.372548997
     Max Dmom(:,3) =  -0.372548997 Min Dmom(:,3) =  -0.372548997
     Max Dener     =  -0.372548997 Min Dener     =  -0.372548997

