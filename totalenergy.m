function tot=totalenergy(kxlist,kylist,t_bond,t,n_bond,U,energyall,wfall,parameters)
N=length(kxlist);
Q=parameters.Q;
NQ=length(Q);
Qindex=parameters.Qindex;
% Qindex=parameters.Qindex;
% Qindexmod=parameters.Qindexmod;
kxbasis=cell(1,NQ);
kybasis=cell(1,NQ);
for i=1:NQ
    kxbasis{i}=kxlist+Q{i}(1);
    kybasis{i}=kylist+Q{i}(2);
end

energylist=real(tb(t_bond,t,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));

energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
occupied=(energyall<=mu);

T=zeros(2*NQ,2*NQ,N);
for k_beta_index=1:N
    T(:,:,k_beta_index)=diag(energylist(k_beta_index,:));
end

Qx=cellfun(@(x)x(1),Q);
Qy=cellfun(@(x)x(2),Q);
[q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);

V1=V(n_bond,U,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y); %V1_{q_alpha,q_delta}

delta_tensor=zeros(NQ,NQ,NQ,NQ); %delta_tensor_{q_alpha,q_beta,q_gamma,q_delta}
for q_alpha_index=1:NQ
    for q_beta_index=1:NQ
        for q_gamma_index=1:NQ
            for q_delta_index=1:NQ
                qindex_alpha=Qindex{q_alpha_index};
                qindex_beta=Qindex{q_beta_index};
                qindex_gamma=Qindex{q_gamma_index};
                qindex_delta=Qindex{q_delta_index};
                deltafunc=qindex_gamma+qindex_delta-qindex_alpha-qindex_beta;
                delta_tensor(q_alpha_index,q_beta_index,q_delta_index,q_gamma_index)=all(deltafunc==round(deltafunc));
            end
        end
    end
end

c=zeros(N,2*NQ,NQ,2);
c(:,:,:,1)=wfall(:,:,1:NQ);  % k_alpha,n,q_alpha,spin_up
c(:,:,:,2)=wfall(:,:,NQ+1:2*NQ); % k_alpha,n,q_alpha,spin_down

epsilonk=zeros(N,NQ,2); %k,q,sigma
epsilonk(:,:,1)=energyall(:,1:NQ);
epsilonk(:,:,2)=energyall(:,NQ+1:2*NQ);

epsilonk_expand=permute(repmat(epsilonk,1,1,1,2*NQ),[1,4,2,3]); %k,n,q,sigma

prod=sum((conj(c).*c).*epsilonk_expand,[3,4]); %k,n
T=ttt(tensor(prod),tensor(occupied),[1,2],[1,2]);

occupied_expand=repmat(occupied,1,1,NQ,2); %k_alpha,n,q_alpha,sigma1
c_tmp=occupied_expand.*conj(c); % k_alpha,n,q_alpha,sigma

ave_alpha=ttt(tensor(c_tmp),tensor(c),[1,2,4],[1,2,4]); % {q_alpha, q_delta}
ave_beta=ave_alpha; %{q_beta,q_gamma}
prod1=V1.*ave_alpha; %{q_alpha,q_delta}
prod2=ttt(tensor(delta_tensor),ave_beta,[3,4],[1,2]);%{q_alpha,q_delta}
H1=ttt(prod1,prod2,[1,2],[1,2])/(N*NQ); 


c_tmp2=occupied_expand.*conj(c); %k_alpha,n,q_alpha,sigma1
c_tmp2_expand=permute(repmat(c_tmp2,1,1,1,1,NQ,2),[1,2,3,5,4,6]); %{k_alpha,n,q_alpha,q_gamma,sigma1,sigma2}
c_expand=permute(repmat(c,1,1,1,1,NQ,2),[1,2,5,3,6,4]); %{k_alpha,n,q_alpha,q_gamma,sigma1,sigma2}
ave2_alpha=squeeze(sum(c_tmp2_expand.*c_expand,2)); %{k_alpha,q_alpha,q_gamma,sigma1,sigma2}
ave2_beta=permute(ave2_alpha,[1,2,3,5,4]); %{k_beta,q_beta,q_delta,sigma1,sigma2}

ave2=permute(ttt(tensor(ave2_alpha),tensor(ave2_beta),[4,5],[5,4]),[1,4,2,5,3,6]); %{k_alpha,k_beta,q_alpha,q_beta,q_gamma,q_delta}

[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
V2=V(n_bond,U,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y); %V2_{k_alpha,k_beta,q_alpha,q_delta}
V2_expand=permute(repmat(V2,1,1,1,1,NQ,NQ),[1,2,3,5,6,4]); %{k_alpha,k_beta,q_alpha,q_beta,q_gamma,q_delta}
prod2=squeeze(sum(V2_expand.*ave2.data,[1,2])); %{q_alpha,q_beta,q_gamma,q_delta}
H2=ttt(tensor(prod2),tensor(delta_tensor),[1,2,3,4],[1,2,3,4])/(N*NQ);
tot=T+real(H1-H2)/2;
end





