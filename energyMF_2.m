function [energyall,wfall]=energyMF_2(kxlist,kylist,t_bond,t,n_bond,U,energyall_o,wfall_o,parameters)
%use tt2
N=length(kxlist);
Q=parameters.Q;
NQ=length(Q);
Qindex=parameters.Qindex;
Qindexmod=parameters.Qindexmod;
kxbasis=cell(1,NQ);
kybasis=cell(1,NQ);
for i=1:NQ
    kxbasis{i}=kxlist+Q{i}(1);
    kybasis{i}=kylist+Q{i}(2);
end

energylist=real(tb(t_bond,t,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));

energyall=zeros(N,2*NQ);
wfall=zeros(N,2*NQ,2*NQ);


Qx=cellfun(@(x)x(1),Q);
Qy=cellfun(@(x)x(2),Q);
[q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);

V1=V(n_bond,U,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}

delta_tensor=zeros(NQ,NQ,NQ,NQ); %delta_tensor_{q_alpha,q_beta,q_gamma,q_delta}
for q_alpha_index=1:NQ
    for q_beta_index=1:NQ
        for q_gamma_index=1:NQ
            for q_delta_index=1:NQ
%                 qindex_alpha=Qindex{q_alpha_index};
%                 qindex_beta=Qindex{q_beta_index};
%                 qindex_gamma=Qindex{q_gamma_index};
%                 qindex_delta=Qindex{q_delta_index};
%                 deltafunc=qindex_gamma+qindex_delta-qindex_alpha-qindex_beta;
%                 delta_tensor(q_alpha_index,q_beta_index,q_delta_index,q_gamma_index)=all(deltafunc==round(deltafunc));
                qindex_alpha=Qindexmod{q_alpha_index};
                qindex_beta=Qindexmod{q_beta_index};
                qindex_gamma=Qindexmod{q_gamma_index};
                qindex_delta=Qindexmod{q_delta_index};
                deltafunc=qindex_gamma+qindex_delta-qindex_alpha-qindex_beta;
                delta_tensor(q_alpha_index,q_beta_index,q_delta_index,q_gamma_index)=all(mod(deltafunc,NQ)==0);
            end
        end
    end
end

ave=average(energyall_o,wfall_o,parameters); %k_alpha,q_alpha,q_delta,sigma1,sigma2
ave1=contract(tensor(ave),4,5); %k_alpha,q_alpha,q_delta
V1_expand=permute(repmat(V1,1,1,NQ,NQ),[1,3,4,2]); %q_alpha,q_beta,q_gamma,q_delta
prod1=V1_expand.*delta_tensor; %q_alpha,q_beta,q_gamma,q_delta
prod2=ttt(ave1,tensor(prod1),[2,3],[1,4]); %k_alpha,q_beta,q_gamma
    prod2=1/2*(prod2+conj(permute(prod2.data,[1,3,2])));
elem1=squeeze(sum(prod2.data,1)); %q_beta,q_gamma
H1=kron(eye(2),elem1);
H1=repmat(H1,1,1,N)/(N*NQ); %{q_beta+sigma,q_gamma+sigma,k_beta}

[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);

V2=V(n_bond,U,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}

ave2=ave; %k_alpha,q_alpha,q_gamma,sigma1,sigma2
prod1=ttt2(V2,ave2,[1],[1],[3],[2]); %q_alpha,k_beta,q_delta,q_gamma,sigma1,sigma2
prod2=ttt2(delta_tensor,prod1,[1,3],[1,4],[4],[3]); %q_delta,q_beta,k_beta,sigma1,sigam2
H2=reshape(permute(prod2,[2,5,1,4,3]),[2*NQ,2*NQ,N])/(N*NQ); %q_beta+sigma2,q_delta+sigma1,k_beta

T=zeros(2*NQ,2*NQ,N);
for k_beta_index=1:N
    T(:,:,k_beta_index)=diag(energylist(k_beta_index,:));
end

H=T+H1-H2;
herr=max(sum(abs(H-conj(permute(H,[2,1,3]))),[1,2]));
assert(herr<1e-12,'hermitian error exceeds');
fprintf("hermitian error: %e\n",herr);
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



