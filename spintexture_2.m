function [spin,gap]=spintexture_2(energyall,wfall,parameters)
energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
occupied=(energyall<=mu);
gap=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2))+1)-energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
inner=parameters.inner;
Q=parameters.Q;
NQ=length(Q);
N=size(energyall,1);
% inner{end+1}=(-1)*parameters.aM1+1*parameters.aM2;
% inner{end+1}=(-2)*parameters.aM1+2*parameters.aM2;
Nai=length(inner);

sigma(:,:,1)=eye(2);
sigma(:,:,2)=[0,1;1,0];
sigma(:,:,3)=[0,-1i;1i,0];
sigma(:,:,4)=[1,0;0,-1];

expqq=zeros(Nai,NQ,NQ); %expqq_{ai,q_alpha,q_delta}
for ai_index=1:Nai
    for q_alpha_index=1:NQ
        for q_delta_index=1:NQ
            q_alpha=Q{q_alpha_index};
            q_delta=Q{q_delta_index};
            expqq(ai_index,q_alpha_index,q_delta_index)=exp(-1i*(inner{ai_index}*(q_alpha-q_delta)'));
        end
    end
end

c=reshape(wfall,N,2*NQ,NQ,2); %k,n,q,sigma
cc=ttt2(conj(c),c,[],[],[1,2],[1,2]); %k,n,q1,sigma1,q2,sigma2
prod1=ttt(tensor(occupied),tensor(cc),[1,2],[1,2]); %q1,sigma1,q2,sigma2
prod2=ttt(prod1,tensor(sigma),[2,4],[1,2]); %q1,q2,alpha
spin=ttt(tensor(expqq),prod2,[2,3],[1,2])/(N*NQ); %ai,alpha

spin=real(spin.data);
end

