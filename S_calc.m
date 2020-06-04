function [S,mu]=spin(energyall,wfall,parameters)
%S(a,b,c), a,b=1(up) or 2(down) c=1,2,3 (0, Q, -Q)

energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/parameters.nu(2));
occupied=(energyall<=mu);
wf_occupied=[wfall{occupied}];

                