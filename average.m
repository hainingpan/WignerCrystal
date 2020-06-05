function ave=average(energyall,wfall,parameters)
%calculate <c_{k1+q_1,sigma_1}^dagger c_{k1+q_2,sigma2}>
energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
occupied=(energyall<=mu);
Q=parameters.Q;
NQ=length(Q);
N=size(energyall,1);

c=reshape(wfall,N,2*NQ,NQ,2); %k,n,q,sigma
    
c_expand=ttt2(conj(c),c,[],[],[1,2],[1,2]);

ave=ttt2(tensor(c_expand),tensor(occupied),[2],[2],[1],[1]);
ave=permute(ave,[1,2,4,3,5]); %k,q1,q2,sigma1,sigma2
end