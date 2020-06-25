function [energyall,wfall]=energyMF(kxlist,kylist,t_bond,t,n_bond,U,energyall_o,wfall_o,parameters)
%Obsolete (wrong), use energyMF_2 instead; for record purpose only

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

energyall=zeros(N,2*NQ);
wfall=zeros(N,2*NQ,2*NQ);
energyall_o_sort=sort(energyall_o(:));
mu=energyall_o_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
occupied=(energyall_o<=mu);

Qx=cellfun(@(x)x(1),Q);
Qy=cellfun(@(x)x(2),Q);
[q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);

V1=V(n_bond,U,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y); %V1_{q_alpha,q_delta}


% V1=zeros(NQ,NQ);    %V1_{q_alpha,q_delta}
% for q_alpha_index=1:NQ 
%     for q_delta_index=1:NQ
%         q_alpha=Q{q_alpha_index};
%         q_delta=Q{q_delta_index};
%         V1(q_alpha_index,q_delta_index)=V(n_bond,U,q_alpha-q_delta);
%     end
% end

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
c(:,:,:,1)=wfall_o(:,:,1:NQ);  % k_alpha,n,q_alpha,sigma1
c(:,:,:,2)=wfall_o(:,:,NQ+1:2*NQ); % k_alpha,n,q_alpha,sigma1

occupied_expand=repmat(occupied,1,1,NQ,2); %k_alpha,n,q_alpha,sigma1
c_tmp=occupied_expand.*conj(c); % k_alpha,n,q_alpha,sigma

ave=ttt(tensor(c_tmp),tensor(c),[1,2,4],[1,2,4])/(N*NQ); % {q_alpha, q_delta}
elem=ttt(tensor(delta_tensor),tensor(V1.*ave.data),[1,4],[1,2]); % {q_beta,q_gamma}
H1=kron(eye(2),elem.data);
H1=repmat(H1,1,1,N); %{q_beta,q_gamma,k_beta}

c_tmp2=occupied_expand.*conj(c);
c_tmp2_expand=permute(repmat(c_tmp2,1,1,1,1,NQ,2),[1,2,3,5,4,6]); %{k_alpha,n,q_alpha,q_gamma,sigma1,sigma2}
c_expand=permute(repmat(c,1,1,1,1,NQ,2),[1,2,5,3,6,4]); %{k_alpha,n,q_alpha,q_gamma,sigma1,sigma2}
ave2=squeeze(sum(c_tmp2_expand.*c_expand,2))/(N*NQ); %{k_alpha,q_alpha,q_gamma,sigma1,sigma2}

delta_tensor_expand=permute(repmat(delta_tensor,1,1,1,1,N,2,2),[5,1,2,3,4,6,7]); %{k_alpha,q_alpha,q_beta,q_gamma,q_delta,sigma1,sigma2}
ave2_expand=permute(repmat(ave2,1,1,1,1,1,NQ,NQ),[1,2,6,3,7,4,5]); %{k_alpha,q_alpha,q_beta,q_gamma,q_delta,sigma1,sigma2}
prod2=squeeze(sum(delta_tensor_expand.*ave2_expand,4)); %{k_alpha,q_alpha,q_beta,q_delta,sigma1,sigma2}

[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);

V2=V(n_bond,U,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y); %V2_{k_alpha,k_beta,q_alpha,q_delta}

elem2=zeros(N,NQ,NQ,2,2); %{k_beta,q_beta,q_delta,sigma1,sigma2}
for q_delta_index=1:NQ 
    elem2(:,:,q_delta_index,:,:)=ttt(tensor(V2(:,:,:,q_delta_index)),tensor(prod2(:,:,:,q_delta_index,:,:)),[1,3],[1,2]);
end


H2=zeros(2*NQ,2*NQ,N);
for k_beta_index=1:N
    H2(:,:,k_beta_index)=[squeeze(elem2(k_beta_index,:,:,1,1)),squeeze(elem2(k_beta_index,:,:,2,1));
                          squeeze(elem2(k_beta_index,:,:,1,2)),squeeze(elem2(k_beta_index,:,:,2,2))];
end


T=zeros(2*NQ,2*NQ,N);
for k_beta_index=1:N
    T(:,:,k_beta_index)=diag(energylist(k_beta_index,:));
end

H=T+H1-H2;
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
%     k_beta=[kxlist(k_beta_index),kylist(k_beta_index)];
%     V_tmp=zeros(N,NQ,NQ); % {k_alpha,q_alpha,q_delta}
%     for k_alpha_index=1:N        
%         for q_alpha_index=1:NQ
%             for q_delta_index=1:NQ
%                 k_alpha=[kxlist(k_alpha_index),kylist(k_alpha_index)];
%                 V_tmp(k_alpha_index,q_alpha_index,q_delta_index)=V(n_bond,U,k_alpha+q_alpha-k_beta-q_delta);
%             end
%         end
%     end 
%     V_tmp_expand=repmat(V_tmp,1,1,1,NQ,NQ); %{k_alpha,q_alpha,q_beta,q_gamma,q_delta}
%     delta_tensor_expand=permute(repmat(delta_tensor,1,1,1,1,N),[5,1,2,3,4]); %{k_alpha,q_alpha,q_beta,q_gamma,q_delta}
%     V2=V_tmp_expand.*delta_tensor_expand;   %{k_alpha,q_alpha,q_beta,q_gamma,q_delta}
%     elem2=ttt(tensor(V2),tensor(ave2),[1,2,4],[1,2,3]); %{q_beta,q_delta,sigma1,sigma2}
%     elem2=elem2.data;
%     H2=[elem2(:,:,1,1),elem2(:,:,2,1);
%         elem2(:,:,1,2),elem2(:,:,2,2)];
    
%     for q_alpha_index=1:NQ
%         for q_beta_index=1:NQ
%             for q_gamma_index=1:NQ
%                 for q_delta_index=1:NQ
%                     qindex_alpha=Qindex{q_alpha_index};
%                     qindex_beta=Qindex{q_beta_index};
%                     qindex_gamma=Qindex{q_gamma_index};
%                     qindex_delta=Qindex{q_delta_index};
%                     deltafunc=qindex_gamma+qindex_delta-qindex_alpha-qindex_beta;
%                     if all(deltafunc==round(deltafunc))                                           
%                         for sigma_1=0:1     %0: spin up; 1: spin down; sigma_1*2*NQ is the shift of index
%                             for sigma_2=0:1
%                                 q_alpha=Q{q_alpha_index};
%                                 q_beta=Q{q_beta_index};
%                                 q_gamma=Q{q_gamma_index};
%                                 q_delta=Q{q_delta_index};
% %                                 bra1=locateindex(q_beta,Qlistmod);
% %                                 ket1=locateindex(q_gamma,Qlistmod);
%                                 V1=V(n_bond,U,q_alpha-q_delta);
%                                 sum1=0;
%                                 sum2=0;
%                                 for k_alpha_index=1:N
%                                     for level=1:2*NQ
%                                         k_alpha=[kxbasis{k_alpha_index},kybasis{k_alpha_index}];
%                                         
%                                         sum1=sum1+conj(wfall_o{k_alpha_index,level}(q_alpha_index+sigma_1*2*NQ))*...
%                                             wfall_o{k_alpha_index,level}(q_delta_index+sigma_1*2*NQ)*occupied(k_alpha_index,level);
% 
%                                         sum2=sum2+V(n_bond,U,k_alpha+q_alpha-k_beta-q_delta)...
%                                             *conj(wfall_o{k_alpha_index,level}(q_alpha_index+sigma_1*2*NQ))*...
%                                             wfall_o{k_alpha_index,level}(q_gamma_index+sigma_2*2*NQ)*occupied(k_alpha_index,level);
%                                     end
%                                 end
%                                 H1(q_beta_index+sigma_2*NQ,q_gamma_index+sigma_2*NQ)=1/N*V1*sum1;
%                                 H2(q_beta_index+sigma_2*NQ,q_delta_index+sigma_1*NQ)=1/N*sum2;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end    
%    H=T+H1-H2; 
%     [vec,val]=eig(H);
%     val=real(diag(val));
%     [val,I]=sort(val);
%     vec=vec(:,I);
%     energyall(k_beta_index,:)=val;
%     for ii=1:2*NQ
%         wfall(k_beta_index,ii,:)=vec(:,ii);
%     end
% end
end



