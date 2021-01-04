%combine V for theta
thetalist=3:0.04:5;
for thetai=1:length(thetalist)    
    load(sprintf('phase4,12_theta(%0.2f,%0.2f,1).mat',thetalist(thetai),thetalist(thetai)));
    V1store{thetai}=V1{1};
    V2store{thetai}=V2{1};
    energyliststore{thetai}=energylist{1};
end
clear V1 V2 energylist
V1=V1store;
V2=V2store;
energylist=energyliststore;
save(sprintf('phase4,12_theta(%0.2f,%0.2f,%d).mat',thetalist(1),thetalist(end),length(thetalist)),'V1','V2','energylist','-v7.3');