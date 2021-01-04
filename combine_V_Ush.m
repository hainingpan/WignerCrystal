%combine V for Ushell
Ushlist=0:1;
for Ushi=1:length(Ushlist)    
    load(sprintf('phase1,1_theta5.00_h1_U%d.mat',Ushlist(Ushi)));
    V1store{Ushi}=V1;
    V2store{Ushi}=V2;
    energyliststore{Ushi}=energylist;
end
clear V1 V2 energylist
V1=V1store;
V2=V2store;
energylist=energyliststore;
save(sprintf('phase1,1_theta5.00_h1_U(%d,%d,%d).mat',Ushlist(1),Ushlist(end),length(Ushlist)),'V1','V2','energylist','-v7.3');