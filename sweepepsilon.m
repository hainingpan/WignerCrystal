function [final,spin,gap,innergap,i]=sweepepsilon(epsilon,parameters)
parameters.V1=parameters.V1/epsilon;
parameters.V2=parameters.V2/epsilon;


[energyall,wfall]=energyMF_init_2(parameters);
    for i=1:10000
       [en(i),ave,V2deltaave]=totalenergy_2(energyall,wfall,parameters);

        if length(en)>1    
            if abs(en(end)-en(end-1))<1e-11
                break
            end
        end
        
        [energyall,wfall]=energyMF_2(ave,V2deltaave,parameters);        
    end
final=en(end);
[spin,gap,innergap]=spintexture(energyall,wfall,parameters);
end
