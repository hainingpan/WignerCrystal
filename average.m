function ave=average(energyall,wfall,parameters)
%calculate <c_{k1+q_1,sigma_1}^dagger c_{k1+q_2,sigma2}>
energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
occupied=(energyall<=mu);
Q=parameters.Q;
NQ=length(Q);
N=size(energyall,1);

c=reshape(wfall,N,2*NQ,NQ,2); %k,n,q,sigma
    
c_expand=tprod(conj(c),c,[],[],[1,2],[1,2]);

ave=tprod((c_expand),(occupied),[2],[2],[1],[1]);
ave=permute(ave,[1,2,4,3,5]); %k,q1,q2,sigma1,sigma2
herr=max(sum(abs(ave-conj(permute(ave,[1,3,2,5,4]))),[2,3]),[],'all');
assert(herr<1e-12,'average spin hermitian error exceeds');
fprintf("average spin hermitian error: %e\n",herr);
ave=1/2*(ave+conj(permute(ave,[1,3,2,5,4])));
end