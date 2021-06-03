%combine V for theta
thetalist=3:0.04:5;
% thetalist=[4,5];
slist=30:10:300;
% slist=[30,60];
V1store={};
V2store={};
energyliststore={};
for thetai=1:length(thetalist)
    tmp=load(sprintf('phase1,1_theta%0.2f_h1_s%d,%d.mat',thetalist(thetai),slist(1),slist(end)));
    for si=1:length(slist)
        V1store{thetai}{si}=tmp.V1{si};
        V2store{thetai}{si}=tmp.V2{si};
        energyliststore{thetai}{si}=tmp.energylist{si};
    end
end
% clear V1 V2 energylist
V1=V1store;
V2=V2store;
energylist=energyliststore;
save(sprintf('phase1,1_theta(%0.2f,%0.2f,%d)_h1_s(%d,%d,%d).mat',thetalist(1),thetalist(end),length(thetalist),slist(1),slist(end),length(slist)),'V1','V2','energylist','-v7.3');