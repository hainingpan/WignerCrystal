function [energyall,wfall]=energyMF_2(ave,V2ave,parameters)
%use ttt2
N=parameters.N;
Q=parameters.Q;
NQ=length(Q);
Ez=parameters.Ez;

energyall=zeros(N,2*NQ);
wfall=zeros(N,2*NQ,2*NQ);

V1=parameters.V1; %V1_{q_alpha,q_delta}
delta_tensor=parameters.delta_tensor; %delta_tensor_{q_alpha,q_beta,q_gamma,q_delta}
% ave=average(energyall_o,wfall_o,parameters); %k_alpha,q_alpha,q_delta,sigma1,sigma2
ave1=contract(tensor(ave),4,5); %k_alpha,q_alpha,q_delta
V1_expand=permute(repmat(V1,1,1,NQ,NQ),[1,3,4,2]); %q_alpha,q_beta,q_gamma,q_delta
prod1=V1_expand.*delta_tensor; %q_alpha,q_beta,q_gamma,q_delta
prod2=ttt(ave1,tensor(prod1),[2,3],[1,4]); %k_alpha,q_beta,q_gamma
prod2=1/2*(prod2+conj(permute(prod2.data,[1,3,2])));
elem1=squeeze(sum(prod2.data,1)); %q_beta,q_gamma
H1=kron(eye(2),elem1);
H1=repmat(H1,1,1,N)/(N*NQ); %{q_beta+sigma,q_gamma+sigma,k_beta}

% V2=parameters.V2; %V2_{k_alpha,k_beta,q_alpha,q_delta}

% ave2=ave; %k_alpha,q_alpha,q_gamma,sigma1,sigma2
% V2ave=ttt2(V2,ave2,[1],[1],[3],[2]); %q_alpha,k_beta,q_delta,q_gamma,sigma1,sigma2
prod2=ttt2(delta_tensor,V2ave,[1,3],[1,4],[4],[3]); %q_delta,q_beta,k_beta,sigma1,sigam2
H2=reshape(permute(prod2,[2,5,1,4,3]),[2*NQ,2*NQ,N])/(N*NQ); %q_beta+sigma2,q_delta+sigma1,k_beta

%q_beta,k_beta,q_delta,sigma1,sigma2
% H2=reshape(permute(V2deltaave,[1,5,3,4,2]),[2*NQ,2*NQ,N])/(N*NQ); %q_beta+sigma2,q_delta+sigma1,k_beta

energylist=parameters.energylist+repmat(Ez*[ones(1,NQ),-ones(1,NQ)],[N,1]);
T=zeros(2*NQ,2*NQ,N);
for k_beta_index=1:N
    T(:,:,k_beta_index)=diag(energylist(k_beta_index,:));
end

H=T+H1-H2;
herr=max(sum(abs(H-conj(permute(H,[2,1,3]))),[1,2]));
assert(herr<1e-12,'hermitian error exceeds');
% fprintf("hermitian error: %e\n",herr);
H=1/2*(H+conj(permute(H,[2,1,3])));
for k_beta_index=1:N
   [vec,val]=eig(H(:,:,k_beta_index));
    val=real(diag(val));
    [val,I]=sort(val);
    vec=vec(:,I);
    energyall(k_beta_index,:)=val;
    for ii=1:2*NQ
        wfall(k_beta_index,ii,:)=vec(:,ii);
    end
end

end



