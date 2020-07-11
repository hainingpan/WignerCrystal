function [tot,ave,V2ave]=totalenergy_2(energyall,wfall,parameters)
N=parameters.N;
Q=parameters.Q;
NQ=length(Q);
Ez=parameters.Ez;

energylist=parameters.energylist;

epsilonk=reshape(energylist+repmat(Ez*[ones(1,NQ),-ones(1,NQ)],[N,1]),[N,NQ,2]); %k_alpha,q_alpha,sigma1

ave=average(energyall,wfall,parameters); %k_alpha,q_alpha,q_delta,sigma1,sigma2
aveT=zeros(N,NQ,2); %k_alpha,q_alpha,sigma1
for q_alpha_index=1:NQ
    for sigma_index=1:2
        aveT(:,q_alpha_index,sigma_index)=ave(:,q_alpha_index,q_alpha_index,sigma_index,sigma_index);
    end
end

T=sum(epsilonk.*aveT,'all'); 

V1=parameters.V1;  %V1_{q_alpha,q_delta}
delta_tensor=parameters.delta_tensor; %delta_tensor_{q_alpha,q_beta,q_gamma,q_delta}
ave1_alpha=squeeze(sum(double(contract(tensor(ave),4,5)),1)); %q_alpha,q_delta
ave1_beta=ave1_alpha; %q_beta,q_gamma
prod1=V1.*ave1_alpha; %q_alpha,q_delta
prod2=ttt(tensor(delta_tensor),tensor(prod1),[1,4],[1,2]); %q_beta,q_gamma
H1=ttt(tensor(ave1_beta),prod2,[1,2],[1,2])/(2*N*NQ);

V2=parameters.V2; %V2_{k_alpha,k_beta,q_alpha,q_delta}
ave2_alpha=ave; %k_alpha,q_alpha,q_gamma,sigma1,sigma2
ave2_beta=ave; %%k_beta,q_beta,q_delta,sigma1,sigma2

% deltaave=ttt2(delta_tensor,ave2_alpha,[3],[3],[1],[2]); %q_alpha,q_beta,q_delta,k_alpha,sigma1,sigma2
% V2deltaave=ttt2(V2,deltaave,[1,3],[4,1],[4],[2]); %q_beta,k_beta,q_delta,sigma1,sigma2
% H2=ttt(tensor(ave2_beta),tensor(V2deltaave),[1,2,3,4,5],[2,1,3,5,4])/(2*N*NQ);
 
V2ave=ttt2(tensor(V2),tensor(ave2_alpha),[1],[1],[3],[2]); %q_alpha,k_beta,q_delta,q_gamma,sigma,sigma'
prod2=ttt2(V2ave,tensor(ave2_beta),[2,5,6],[1,5,4],[3],[3]); %q_delta,q_alpha,q_gamma,q_beta
H2=ttt(tensor(delta_tensor),tensor(prod2),[1,2,3,4],[2,4,3,1])/(2*N*NQ);

tot=real(T+H1-H2);
tot=tot/(N*NQ);
end




