function [final,spin,gap,innergap,i,ch]=sweepepsilon(epsilon,kxlist,kylist,parameters)
parameters.V1=parameters.V1/epsilon;
parameters.V2=parameters.V2/epsilon;

if parameters.nu==[4,8] | parameters.nu==[5,10]
    [ave,V2ave]=average_kagome(parameters.phi,parameters.s1,kxlist,kylist,parameters);
    [energyall,wfall]=energyMF_2(ave,V2ave,parameters);
elseif parameters.nu==[4,12] | parameters.nu==[4,6]
    [ave,V2ave]=average_honeycomb(kxlist,kylist,parameters);
    [energyall,wfall]=energyMF_2(ave,V2ave,parameters);
else
    [energyall,wfall]=energyMF_init_2(parameters);
end

    for i=1:10000
       [en(i),ave,V2deltaave]=totalenergy_2(energyall,wfall,parameters);

        if length(en)>1    
            if abs(en(end)-en(end-1))<1e-15
                break
            end
        end
        
        [energyall,wfall]=energyMF_2(ave,V2deltaave,parameters);        
    end
final=en(end);
[spin,gap,innergap]=spintexture(energyall,wfall,parameters);
ch=chern(wfall,parameters);
end
