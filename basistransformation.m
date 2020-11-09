function [energyall,wfall]=basistransformation(kxlist,kylist,parameters)
%% buggy
%basis transformation C_alpha->C_Q

Nk=length(kxlist);
Q=parameters.Q;
NQ=length(parameters.Q);
ave=zeros(NQ,NQ,Nk); %alpha',alpha,k
wf=zeros(Nk,NQ,NQ);
wfall=zeros(Nk,2*NQ,2*NQ);
energyall=zeros(Nk,2*NQ)+100;

r={[0,0];[-parameters.aM1];[parameters.aM2-parameters.aM1];[parameters.aM2-2*parameters.aM1]};
% r={[0,0];[-parameters.aM1];[parameters.aM2-parameters.aM1]};
expalphaQ=zeros(NQ); %q1, alpha'
for Qindex=1:NQ
    for alphaindex=1:NQ
        expalphaQ(Qindex,alphaindex)=exp(-1i*r{alphaindex}*Q{Qindex}')/sqrt(NQ);
    end
end

for kindex=1:Nk
    kx=kxlist(kindex);
    ky=kylist(kindex);
    [ave(:,:,kindex),wf(kindex,:,:),energyall(kindex,1:NQ)]=kagome_h(pi/6,1,[kx,ky],parameters);
%     honeycomb_h(pi/2,[kx,ky],parameters);
    wfall(kindex,1:NQ,1:NQ)=expalphaQ*squeeze(wf(kindex,:,:));
end


end
        

