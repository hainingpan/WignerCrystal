function ave=average_init(parameters)
%calculate <c_{k1+q_1,sigma_1}^dagger c_{k1+q_2,sigma2}>
%=1/NQ*sum_ai <c_{ai,sigma1}^dagger c_{ai,sigma2} exp(1i*ai*(q1-q2))>
Q=parameters.Q;
NQ=length(Q);
inner=parameters.inner;
sigma_x=[0,1;1,0];
sigma_y=[0,-1i;1i,0];
sigma_z=[1,0;0,-1];
spin0=parameters.spin0;
Nai=length(spin0);


cc=zeros(2,2,Nai); %cc_{sigma1,sigma2,ai}
for ai_index=1:Nai
    avespin=conj(1/2*(spin0{ai_index}(1)*sigma_x+spin0{ai_index}(2)*sigma_y+spin0{ai_index}(3)*sigma_z+eye(2)));
    cc(:,:,ai_index)=avespin;
end

expqq=zeros(Nai,NQ,NQ); %expqq_{ai,q_1,q_2}
for ai_index=1:Nai
    for q_alpha_index=1:NQ
        for q_delta_index=1:NQ
            q_alpha=Q{q_alpha_index};
            q_delta=Q{q_delta_index};
            expqq(ai_index,q_alpha_index,q_delta_index)=exp(1i*(inner{ai_index}*(q_alpha-q_delta)'));
        end
    end
end

ave=ttt(tensor(expqq),tensor(cc,[2,2,Nai]),[1],[3]); %ave_{q_1,q_2,sigma1,sigma2}
ave=ave.data/NQ;

herr=max(sum(abs(ave-conj(permute(ave,[2,1,4,3]))),[2,1]),[],'all');
assert(herr<1e-12,'init average spin hermitian error exceeds');
% fprintf("init average spin hermitian error: %e\n",herr);
ave=1/2*(ave+conj(permute(ave,[2,1,4,3])));

end
