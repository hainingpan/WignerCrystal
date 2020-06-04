function spin=spintexture(energyall,wfall,parameters)
energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
occupied=(energyall<=mu);
inner=parameters.inner;
Q=parameters.Q;
NQ=length(Q);
N=size(energyall,1);
Nai=length(parameters.spin0);

sigma(:,:,1)=eye(2);
sigma(:,:,2)=[0,1;1,0];
sigma(:,:,3)=[0,1i;-1i,0];
sigma(:,:,4)=[1,0;0,-1];

c=zeros(N,2*NQ,NQ,2);
c(:,:,:,1)=wfall(:,:,1:NQ);
c(:,:,:,2)=wfall(:,:,NQ+1:2*NQ);

occupied_expand=repmat(occupied,1,1,NQ,2);
c_tmp=occupied_expand.*conj(c);

c2=ttt(tensor(c),tensor(sigma),[4],[2]);

ave=ttt(tensor(c_tmp),c2,[1,2,4],[1,2,4])/(N*NQ);

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

spin=ttt(tensor(expqq),ave,[2,3],[1,2]);
spin=spin.data;


% for ailist=1:length(inner)
%     for sigma_index=1:4
%         sum=0;
%         for k1_index=1:N
%             for q1_index=1:NQ
%                 for q2_index=1:NQ
%                     q1=Q{q1_index};
%                     q2=Q{q2_index};
%                     for level=1:2*NQ
%                         bra=[wfall{k1_index,level}(q1_index);wfall{k1_index,level}(q1_index+NQ)];
%                         ket=[wfall{k1_index,level}(q2_index);wfall{k1_index,level}(q2_index+NQ)];
%                         sum=sum+bra'*sigma{sigma_index}*ket*exp(-1i*(q1-q2)*inner{ailist}')*occupied(k1_index,level);
%                     end
%                 end
%             end            
%         end
%         sum=sum/(N*NQ);
%         spin(ailist,sigma_index)=sum;
%     end
end

