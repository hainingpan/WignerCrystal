function re=phasediagram(theta,Vz,epsilonlist)

parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'nu',[1,1],'d',inf,'Vz',Vz);
tshell=3;
Ushell=0;
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

t_bond=[neighborlist{1:tshell+1}];
U_bond=[neighborlist{1:Ushell+1}];
hp=1;
tlist=hp*[t{1:tshell+1}];
re=zeros(1,length(epsilonlist));

for epi=1:length(epsilonlist)

    epsilon=epsilonlist(epi);
    Ulist=real([U{1:Ushell+1}])/epsilon;

    clear spinsav en
    [energyall,wfall]=energyMF_init_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,parameters);
    for i=1:100
    [~,gap]=spintexture(energyall,wfall,parameters);
    en(i)=totalenergy_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,energyall,wfall,parameters);
    if length(en)>1    
        if abs(en(end)-en(end-1))<1e-5
            break
        end
    end
    [energyall,wfall]=energyMF_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,energyall,wfall,parameters);
    end
    re(epi)=gap;
end