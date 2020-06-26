function sum=V(nlist,Ulist,kx,ky,parameters)
sum=0;
parfor i=1:length(Ulist)
    n=nlist{i}(1)*parameters.aM1+nlist{i}(2)*parameters.aM2;
    sum=sum+Ulist(i)*exp(-1i*(kx.*n(1)+ky.*n(2)));
end
sum=real(sum);
end
