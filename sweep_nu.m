function sweep_nu(nu,d,epsilon)
% nu=[1,3];
% d=10;
% epsilon=10;
parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',nu,'d',d);
tshell=3;
Ushell=110;
[t,neighborlist]=t_calc_func(tshell,parameters);
U=U_calc_func_2(Ushell,parameters);


n=15;
counter=1;
kxlist=zeros(1,n^2);
kylist=zeros(1,n^2);
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

t_bond=[neighborlist{1:tshell+1}];
U_bond=[neighborlist{1:Ushell+1}];
hp=1;
tlist=-hp*[t{1:tshell+1}];
Ulist=real([U{1:Ushell+1}])/epsilon;

clear spinsav en
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
[spin,gap]=spintexture(energyall,wfall,parameters);

save(sprintf('nu%d,%d_d%d_ep%d_U%d.mat',parameters.nu(1),parameters.nu(2),d,epsilon,Ushell),'final','gap','spin','nu')
