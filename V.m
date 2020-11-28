function sum=V(nlist,Ulist,kx,ky,parameters)
%assume the first element is U0
% sum=Ulist(1)*ones(size(kx));
sum=Ulist(1);
% sum=0;
for i=2:length(Ulist)
    n=nlist{i}(1)*parameters.aM1+nlist{i}(2)*parameters.aM2;
%     if (n(1))==0 & (n(2))==0
%         sum=sum+Ulist(i)*exp(-1i*(kx.*0+ky.*0));
%     else
    if n(2)>1e-12 | abs(n(2))<1e-15 & n(1)>-1e-12
        sum=sum+Ulist(i)*2*cos(kx.*n(1)+ky.*n(2));
    end
%      sum=sum+Ulist(i)*exp(-1i*(kx.*n(1)+ky.*n(2)));
end
sum=real(sum);
end
