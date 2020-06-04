function sum=V(nlist,Ulist,kx,ky)
sum=0;
for i=1:length(Ulist)
    sum=sum+Ulist(i)*exp(-1i*(kx.*nlist{i}(1)+ky.*nlist{i}(2)));
end
