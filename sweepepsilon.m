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
