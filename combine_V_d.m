%combine V for theta
% thetalist=3:0.04:5;
thetalist=[4,5];
% dlist=30:10:300;
dlist=[30,60];
V1store={};
V2store={};
energyliststore={};
for thetai=1:length(thetalist)
    tmp=load(sprintf('phase1,1_theta%0.2f_h1_d%d,%d.mat',thetalist(thetai),dlist(1),dlist(end)));
    for di=1:length(dlist)
        V1store{thetai}{di}=tmp.V1{di};
        V2store{thetai}{di}=tmp.V2{di};
        energyliststore{thetai}{di}=tmp.energylist{di};
    end
end
% clear V1 V2 energylist
V1=V1store;
V2=V2store;
energylist=energyliststore;
save(sprintf('phase1,1_theta(%0.2f,%0.2f,%d)_h1_d(%d,%d,%d).mat',thetalist(1),thetalist(end),length(thetalist),dlist(1),dlist(end),length(dlist)),'V1','V2','energylist','-v7.3');