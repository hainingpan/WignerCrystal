function ave=average(energyall,wfall,parameters)
%calculate <c_{k1+q_1,sigma_1}^dagger c_{k1+q_2,sigma2}>
energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
if parameters.T==0
    occupied=(energyall<=mu); %k,n
else
    occupied=1./(1+exp((energyall-mu)/(parameters.T)));
end
Q=parameters.Q;
NQ=length(Q);
N=size(energyall,1);

c=reshape(wfall,N,2*NQ,NQ,2); %k,n,q,sigma
    
% c_expand=ttt2(conj(c),c,[],[],[1,2],[1,2]);
prod1=repmat(occupied,[1,1,NQ,2]).*c;
prod1=permute(repmat(prod1,[1,1,1,1,NQ,2]),[1,2,5,3,6,4]);
cconj=permute(repmat(conj(c),[1,1,1,1,NQ,2]),[1,2,3,5,4,6]);
ave=squeeze(sum(prod1.*cconj,2));
% ave=ttt2(conj(c),prod1,[2],[2],[1],[1]);

% ave=ttt2(tensor(c_expand),tensor(occupied),[2],[2],[1],[1]);
% ave=permute(ave,[1,2,4,3,5]); %k,q1,q2,sigma1,sigma2
herr=max(sum(abs(ave-conj(permute(ave,[1,3,2,5,4]))),[2,3]),[],'all');
assert(herr<1e-12,'average spin hermitian error exceeds');
% fprintf("average spin hermitian error: %e\n",herr);
ave=1/2*(ave+conj(permute(ave,[1,3,2,5,4])));
end