function [final,spin,gap,innergap]=sweepepsilon(epsilon,parameters)
parameters.V1=parameters.V1/epsilon;
parameters.V2=parameters.V2/epsilon;


[energyall,wfall]=energyMF_init_2(parameters);
    for i=1:500
        en(i)=totalenergy_2(energyall,wfall,parameters);
%         [spin,gap]=spintexture(energyall,wfall,parameters);
%         fprintf("%d: gap:%0.8f meV E:%f meV\n",i,1000*gap,1000*en(end));
%         disp([spin,angle(spin(:,2)+spin(:,3)*1i)*180/pi,angle(spin(:,4)+sqrt(spin(:,3).^2+spin(:,2).^2)*1i)*180/pi]); 
        [energyall,wfall]=energyMF_2(energyall,wfall,parameters);
        if length(en)>1    
            if abs(en(end)-en(end-1))<1e-11
                break
            end
        end
    end
final=en(end);
[spin,gap,innergap]=spintexture(energyall,wfall,parameters);
end
