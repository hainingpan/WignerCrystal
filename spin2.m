function spin=spin2(ave,parameters)
%% buggy
Q=parameters.Q;
NQ=length(Q);
inner=parameters.inner;
N=size(ave,1);
Nai=length(inner);

ave=squeeze(sum(ave,1));
% ave=permute(reshape(permute(ave,[3,4,1,2]),[4,NQ,NQ]),[2,3,1]);
ave2=zeros(NQ,NQ,4);
for i=1:NQ
    for j=1:NQ
        zz=ave(i,j,:,:);
        ave2(i,j,:)=zz(:);
    end
end
ave=ave2;

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

spin=ttt(tensor(expqq),tensor(ave),[2,3],[1,2])/(N*NQ);
spin=real(spin.data);

end