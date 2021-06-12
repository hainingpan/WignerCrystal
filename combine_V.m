function combine_V(nu,Ushell)
%combine V for theta
thetalist=3:0.04:5;
for thetai=1:length(thetalist)        
    load(sprintf('phase%d,%d_theta(%0.2f,%0.2f,1)_h(1)_U%d.mat',nu(1),nu(2),thetalist(thetai),thetalist(thetai),Ushell));
    V1store{thetai}=V1{1};
    V2store{thetai}=V2{1};
    energyliststore{thetai}=energylist{1};
end
clear V1 V2 energylist
V1=V1store;
V2=V2store;
energylist=energyliststore;
save(sprintf('phase%d,%d_theta(%0.2f,%0.2f,%d)_h(1)_U%d.mat',nu(1),nu(2),thetalist(1),thetalist(end),length(thetalist),Ushell),'V1','V2','energylist','-v7.3');