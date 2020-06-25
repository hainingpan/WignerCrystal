function tot=totalenergy_2(kxlist,kylist,t_bond,t,n_bond,U,energyall,wfall,parameters)
N=length(kxlist);
Q=parameters.Q;
NQ=length(Q);
% Qindex=parameters.Qindex;
% Qindexmod=parameters.Qindexmod;
kxbasis=cell(1,NQ);
kybasis=cell(1,NQ);
for i=1:NQ
    kxbasis{i}=kxlist+Q{i}(1);
    kybasis{i}=kylist+Q{i}(2);
end

energylist=real(tb(t_bond,t,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters)); %k,q+sigma
epsilonk=reshape(energylist,[N,NQ,2]);

ave=average(energyall,wfall,parameters); %k_alpha,q_alpha,q_delta,sigma1,sigma2
aveT=zeros(N,NQ,2);
for q_alpha_index=1:NQ
    for sigma_index=1:2
        aveT(:,q_alpha_index,sigma_index)=ave(:,q_alpha_index,q_alpha_index,sigma_index,sigma_index);
    end
end

T=sum(epsilonk.*aveT,'all'); 

Qx=cellfun(@(x)x(1),Q);
Qy=cellfun(@(x)x(2),Q);
[q_alpha_x,q_delta_x]=ndgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=ndgrid(Qy,Qy);

V1=V(n_bond,U,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}

delta_tensor=parameters.delta_tensor; %delta_tensor_{q_alpha,q_beta,q_gamma,q_delta}

ave1_alpha=squeeze(sum(double(contract(tensor(ave),4,5)),1)); %q_alpha,q_delta
ave1_beta=ave1_alpha; %q_beta,q_gamma
prod1=V1.*ave1_alpha; %q_alpha,q_delta
prod2=ttt(tensor(delta_tensor),tensor(prod1),[1,4],[1,2]); %q_beta,q_gamma
H1=ttt(tensor(ave1_beta),prod2,[1,2],[1,2])/(2*N*NQ);

[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);

V2=V(n_bond,U,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}

% ave2_alpha=ave; %k_alpha,q_alpha,q_gamma,sigma1,sigma2
% ave2_beta=ave; %%k_beta,q_beta,q_delta,sigma2,sigma1
% prod1=ttt(tensor(ave2_alpha),tensor(ave2_beta),[4,5],[5,4]); %k_alpha,q_alpha,q_gamma,k_beta,q_beta,q_delta
% prod2=ttt2(tensor(V2),prod1,[1,2],[1,4],[3,4],[2,6]); %q_alpha,q_delta,q_gamma,q_beta
% H2=ttt(tensor(delta_tensor),tensor(prod2),[1,2,3,4],[1,4,3,2])/(2*N*NQ);

ave2_alpha=ave; %k_alpha,q_alpha,q_gamma,sigma1,sigma2
ave2_beta=permute(ave,[1,2,3,5,4]); %%k_beta,q_beta,q_delta,sigma2,sigma1
prod1=ttt2(tensor(V2),tensor(ave2_alpha),[1],[1],[3],[2]); %q_alpha,k_beta,q_delta,q_gamma,sigma,sigma'
prod2=ttt2(prod1,tensor(ave2_beta),[2,5,6],[1,5,4],[3],[3]); %q_delta,q_alpha,q_gamma,q_beta
H2=ttt(tensor(delta_tensor),tensor(prod2),[1,2,3,4],[2,4,3,1])/(2*N*NQ);

tot=real(T+H1-H2);
tot=tot/(N*NQ);
end





