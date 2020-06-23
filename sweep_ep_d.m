%Sweep for Wigner Crystal as a function epsilon and d
dlist=linspace(1,10,20);
Nd=length(dlist);
epsilonlist=linspace(1,60,60);
nu=[1,3];
Nep=length(epsilonlist);
final=zeros(Nd,Nep);
gap=zeros(Nd,Nep);
innergap=zeros(Nd,Nep);

n=15;
kxlist=zeros(1,n^2);
kylist=zeros(1,n^2);
parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',nu);

counter=1;
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

parfor di=1:Nd
parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'d',dlist(di),'nu',nu);
tshell=3;
Ushell=35;
[t,neighborlist]=t_calc_func(tshell,parameters);
U=U_calc_func_2(Ushell,parameters);
sweepfunc=@(x,y) sweepepsilon(tshell,x,y,neighborlist,t,U,kxlist,kylist,parameters);
    for epi=1:Nep
        [final(di,epi),spin(:,:,di,epi),gap(di,epi),innergap(di,epi)]=sweepfunc(Ushell,epsilonlist(epi));
    end
end
nu=parameters.nu;
save('sweepWC.mat','nu','final','spin','epsilonlist','gap','innergap','dlist');


function [final,spin,gap,innergap]=sweepepsilon(tshell2,Ushell2,epsilon,neighborlist,t,U,kxlist,kylist,parameters)
t_bond=[neighborlist{1:tshell2+1}];
U_bond=[neighborlist{1:Ushell2+1}];
hp=1;
tlist=-hp*[t{1:tshell2+1}];
Ulist=real([U{1:Ushell2+1}])/epsilon;

[energyall,wfall]=energyMF_init_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,parameters);
    for i=1:200
        en(i)=totalenergy_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,energyall,wfall,parameters);
        [energyall,wfall]=energyMF_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,energyall,wfall,parameters);
        if length(en)>1    
            if abs(en(end)-en(end-1))<1e-5
                break
            end
        end
    end
final=en(end);
[spin,gap,innergap]=spintexture(energyall,wfall,parameters);
end