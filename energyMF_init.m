function [energyall,wfall]=energyMF_init(kxlist,kylist,t_bond,t,n_bond,U,parameters)
N=length(kxlist);
Q=parameters.Q;
Qindex=parameters.Qindex;
NQ=length(Q);
inner=parameters.inner;
spin0=parameters.spin0;
Nai=length(spin0);
kxbasis=cell(1,NQ);
kybasis=cell(1,NQ);
for i=1:NQ
    kxbasis{i}=kxlist+Q{i}(1);
    kybasis{i}=kylist+Q{i}(2);
end
sigma_x=[0,1;1,0];
sigma_y=[0,-1i;1i,0];
sigma_z=[1,0;0,-1];
energylist=real(tb(t_bond,t,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));

energyall=zeros(N,2*NQ);
wfall=zeros(N,2*NQ,2*NQ);


Qx=cellfun(@(x)x(1),Q);
Qy=cellfun(@(x)x(2),Q);
[q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);

V1=V(n_bond,U,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}

% V1=zeros(NQ,NQ);    %V1_{q_alpha,q_delta}
% for q_alpha_index=1:NQ 
%     for q_delta_index=1:NQ
%         q_alpha=Q{q_alpha_index};
%         q_delta=Q{q_delta_index};
%         V1(q_alpha_index,q_delta_index)=V(n_bond,U,q_alpha(1)-q_delta(1),q_alpha(2)-q_delta(2));
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

cc=zeros(2,Nai); %cc_{sigma1,ai}
for ai_index=1:Nai
    avespin=conj(1/2*(spin0{ai_index}(1)*sigma_x+spin0{ai_index}(2)*sigma_y+spin0{ai_index}(3)*sigma_z+eye(2)));
    cc(:,ai_index)=diag(avespin);
end
cc=(cc)/NQ;

expqq=zeros(Nai,NQ,NQ); %expqq_{ai,q_alpha,q_delta}
for ai_index=1:Nai
    for q_alpha_index=1:NQ
        for q_delta_index=1:NQ
            q_alpha=Q{q_alpha_index};
            q_delta=Q{q_delta_index};
            expqq(ai_index,q_alpha_index,q_delta_index)=exp(1i*(inner{ai_index}*(q_alpha-q_delta)'));
        end
    end
end
ave=ttt(tensor(expqq),tensor(cc),[1],[2]); %ave_{q_alpha,q_delta,sigma}
ave_tmp=sum(ave.data,3)/NQ; %ave_tmp_{q_alpha,q_delta}
% V1_tmp=repmat(V1,1,1,2); %V1_tmp_{q_alpha,q_delta,sigma}
% V1_tmp=ttt(V1,tensor([1,1]),[3],[1]); 
prod=ttt(tensor(delta_tensor),tensor(V1.*ave_tmp),[1,4],[1,2]);
elem=prod.data/NQ;
H1=kron(eye(2),elem); %H1_{q_beta,q_gamma,k_beta}
H1=repmat(H1,1,1,N);

cc2=zeros(2,2,Nai); %cc_{sigma1,ai}
for ai_index=1:Nai
    avespin=conj(1/2*(spin0{ai_index}(1)*sigma_x+spin0{ai_index}(2)*sigma_y+spin0{ai_index}(3)*sigma_z+eye(2)));
    cc2(:,:,ai_index)=avespin;
end
cc2=tensor(cc2)/NQ;
ave2_tmp=ttt(tensor(expqq),cc2,[1],[3]); %ave_tmp_{q_alpha,q_gamma,sigma1,sigma2}
ave2_tmp_expand=permute(repmat(ave2_tmp.data,1,1,1,1,NQ,NQ),[1,5,6,2,3,4]); %{q_alpha,q_beta,q_gamma,q_delta,sigma1,sigma2}
delta_tensor_expand=repmat(delta_tensor,1,1,1,1,2,2); %{q_alpha,q_beta,q_gamma,q_delta,sigma1,sigma2}
prod2=ave2_tmp_expand.*delta_tensor_expand; %{q_alpha,q_beta,q_gamma,q_delta,sigma1,sigma2}

[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);

V2=V(n_bond,U,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}

V2_expand=permute(repmat(squeeze(sum(V2,1)),1,1,1,NQ,2,2),[1,2,4,3,5,6]); %{k_beta,q_alpha,q_beta,q_delta,sigma1,sigma2}
prod2_expand=permute(repmat(squeeze(sum(prod2,3)),1,1,1,1,1,N),[6,1,2,3,4,5]); %{k_beta,q_alpha,q_beta,q_delta,sigma1,sigma2}

elem2=squeeze(sum(V2_expand.*prod2_expand,2))/(N*NQ);
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

   
    
    
    
%     for k_alpha_index=1:N        
%         for q_alpha_index=1:NQ
%             for q_delta_index=1:NQ
%                 k_alpha=[kxlist(k_alpha_index),kylist(k_alpha_index)];
%                 V_tmp(k_alpha_index,q_alpha_index,q_delta_index)=V(n_bond,U,k_alpha+q_alpha-k_beta-q_delta);
%             end
%         end
%     end 
%     
%     
%     V2=squeeze(sum(V_tmp,1)/(N*NQ)); %V2_{q_alpha,q_delta}
%     V2_tmp=permute(repmat(V2,1,1,NQ,NQ,2,2),[1,4,2,3,5,6]);
%     elem2=squeeze(sum(V2_tmp.*prod2,[1,3]));
%     H2=[elem2(:,:,1,1),elem2(:,:,2,1);
%         elem2(:,:,1,2),elem2(:,:,2,2)];
%     
%     
% %     for q_alpha_index=1:NQ
% %         for q_beta_index=1:NQ
% %             for q_gamma_index=1:NQ
% %                 for q_delta_index=1:NQ
% %                     qindex_alpha=Qindex{q_alpha_index};
% %                     qindex_beta=Qindex{q_beta_index};
% %                     qindex_gamma=Qindex{q_gamma_index};
% %                     qindex_delta=Qindex{q_delta_index};
% %                     deltafunc=qindex_gamma+qindex_delta-qindex_alpha-qindex_beta;
% %                     if all(deltafunc==round(deltafunc))                                           
% %                         for sigma_1=0:1     %0: spin up; 1: spin down; sigma_1*2*NQ is the shift of index
% %                             for sigma_2=0:1
% %                                 q_alpha=Q{q_alpha_index};
% %                                 q_beta=Q{q_beta_index};
% %                                 q_gamma=Q{q_gamma_index};
% %                                 q_delta=Q{q_delta_index};
% %                                 bra1=locateindex(q_beta,Qindexmod);
% %                                 ket1=locateindex(q_gamma,Qindexmod);
% %                                 V1=V(n_bond,U,q_alpha-q_delta);
% %                                 sum1=0;
% %                                 sum2=0;
% %                                 
% %                                 for k_alpha_index=1:N
% %                                     k_alpha=[kxlist(k_alpha_index),kylist(k_alpha_index)];
% %                                     for ailist=1:length(spin0)
% %                                         nsigma=1/2*(spin0{ailist}(1)*sigma_x+spin0{ailist}(2)*sigma_y+spin0{ailist}(3)*sigma_z+eye(2));
% %                                         sum1=sum1+nsigma(sigma_1+1,sigma_1+1)*exp(1i*inner{ailist}*(q_alpha-q_delta)');
% %                                         sum2=sum2+V(n_bond,U,k_alpha+q_alpha-k_beta-q_delta)...
% %                                             *nsigma(sigma_1+1,sigma_2+1)*exp(1i*inner{ailist}*(q_alpha-q_gamma)');                                 
% %                                     end
% %                                     sum1=sum1/NQ;
% %                                     sum2=sum2/NQ;
% %                                 end
% %                                 H1(q_beta_index+sigma_2*NQ,q_gamma_index+sigma_2*NQ)=1/N*V1*sum1;
% %                                 H2(q_beta_index+sigma_2*NQ,q_delta_index+sigma_1*NQ)=1/N*sum2;
% %                             end
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% %     end    
%     H=T+H1-H2; 
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


    


