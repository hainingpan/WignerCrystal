function tot=totalenergy_2(kxlist,kylist,t_bond,t,n_bond,U,energyall,wfall,parameters)
N=length(kxlist);
Q=parameters.Q;
NQ=length(Q);
Qindex=parameters.Qindex;
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
[q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);

V1=V(n_bond,U,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}

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

ave1_alpha=squeeze(sum(double(contract(tensor(ave),4,5)),1)); %q_alpha,q_delta
ave1_beta=ave1_alpha; %q_beta,q_gamma
prod1=V1.*ave1_alpha; %q_alpha,q_delta
prod2=ttt(tensor(delta_tensor),tensor(prod1),[1,4],[1,2]); %q_beta,q_gamma
H1=ttt(tensor(ave1_beta),prod2,[1,2],[1,2])/(2*N*NQ);

[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);

V2=V(n_bond,U,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}

ave2_alpha=ave;
ave2_beta=ave;
prod1=ttt(tensor(ave2_alpha),tensor(ave2_beta),[4,5],[5,4]); %k_alpha,q_alpha,q_gamma,k_beta,q_beta,q_delta
prod2=ttt2(tensor(V2),prod1,[1,2],[1,4],[3,4],[2,6]); %q_alpha,q_delta,q_gamma,q_beta
H2=ttt(tensor(delta_tensor),tensor(prod2),[1,2,3,4],[1,4,3,2])/(2*N*NQ);

tot=real(T+H1-H2);
tot=tot/(N*NQ);
end





