function combine_V_vz(nu,Ushell)
%combine V for theta
thetalist=3:0.04:5;
% thetalist=[4,5];
vzlist=1:200;
V1store={};
V2store={};
energyliststore={};
for thetai=1:length(thetalist)
    tmp=load(sprintf('phase%d,%d_theta%0.2f_h1_vz%d,%d_U%d.mat',nu(1),nu(2),thetalist(thetai),vzlist(1),vzlist(end),Ushell));
    for si=1:length(vzlist)
        V1store{thetai}{si}=tmp.V1{si};
        V2store{thetai}{si}=tmp.V2{si};
        energyliststore{thetai}{si}=tmp.energylist{si};
    end
end
% clear V1 V2 energylist
V1=V1store;
V2=V2store;
energylist=energyliststore;
save(sprintf('phase%d,%d_theta(%0.2f,%0.2f,%d)_h1_vz(%d,%d,%d)_U%d.mat',nu(1),nu(2),thetalist(1),thetalist(end),length(thetalist),vzlist(1),vzlist(end),length(vzlist),Ushell),'V1','V2','energylist','-v7.3');