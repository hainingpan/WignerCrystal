function [final,spin,gap,innergap,i,ch]=sweepd(epsilon,kxlist,kylist,parameters)
 parameters.V1=parameters.V1/epsilon;
 parameters.V2=parameters.V2/epsilon;

if ismember(parameters.nu,[[4,8];[5,10];[2,2];[4,4];[7,28];[8,32];[21,28];[24,32]],'row')
    [ave,V2ave]=average_kagome(parameters.phi,parameters.s1,kxlist,kylist,parameters);
    [energyall,wfall]=energyMF_2(ave,V2ave,parameters);
else
    if parameters.nu==[4,12]
        [ave,V2ave]=average_honeycomb(kxlist,kylist,parameters);
        [energyall,wfall]=energyMF_2(ave,V2ave,parameters);
    else
        if ismember(parameters.nu,[[4,6];[6,6]],'row')
            [ave,V2ave]=average_honeycomb2(kxlist,kylist,parameters);
            [energyall,wfall]=energyMF_2(ave,V2ave,parameters);        
        else
            [energyall,wfall]=energyMF_init_2(parameters);
        end
    end
end

for i=1:1000
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
