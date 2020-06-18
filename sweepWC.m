%Sweep for Wigner Crystal

parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',[2,3]);
tshell=3;
Ushell=23;
[t,neighborlist]=t_calc_func(tshell,parameters);
U=U_calc_func(Ushell,parameters);

n=15;
counter=1;
clear kxlist kylist
for xindex=1:n
    for yindex=1:n
        ux=(2*xindex-n-1)/(2*n);
        uy=(2*yindex-n-1)/(2*n);
        klist=ux*parameters.bm1+uy*parameters.bm2;
        kxlist(counter)=klist(1);
        kylist(counter)=klist(2);
        counter=counter+1;
    end
end
kxlist=kxlist';
kylist=kylist';
sweepfunc=@(x,y) sweepepsilon(tshell,x,y,neighborlist,t,U,kxlist,kylist,parameters);
epsilonlist=linspace(1,50,4);
Nep=length(epsilonlist);
final=zeros(1,Nep);
gap=zeros(1,Nep);

parfor i=1:Nep
    disp(i)
[final(i),spin(:,:,i),gap(i)]=sweepfunc(Ushell,epsilonlist(i));
end
save('sweepWC.mat','parameters','final','spin','U','t','epsilonlist','gap');

function [final,spin,gap]=sweepepsilon(tshell2,Ushell2,epsilon,neighborlist,t,U,kxlist,kylist,parameters)
t_bond=[neighborlist{1:tshell2+1}];
U_bond=[neighborlist{1:Ushell2+1}];
hp=1;
tlist=-hp*[t{1:tshell2+1}];
Ulist=real([U{1:Ushell2+1}])/epsilon;

[energyall,wfall]=energyMF_init_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,parameters);
for i=1:200
[spin,gap]=spintexture(energyall,wfall,parameters);
en(i)=totalenergy_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,energyall,wfall,parameters);
[energyall,wfall]=energyMF_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,energyall,wfall,parameters);
if length(en)>1    
    if abs(en(end)-en(end-1))<1e-8
        break
    end
end
end
final=en(end);
end