function re=Vint(rn,epsilon,kx,ky,parameters)
alpha=0.00729735; %e^2/(4*pi*epsilon0); c.f. Quicks Notes on oneNotes
prefactor=alpha/epsilon*2*pi/(sqrt(3)/2*parameters.aM^2);
p=kx-ky;
p=sqrt(p.*p);
plist=p(:);
re=zeros(length(plist),1);
for qindex=1:length(plist)
    q=plist(qindex);
    func=@(r) (1-r./sqrt(r.^2+parameters.d^2)).*besselj(0,q*r);
    re(qindex)=integral(func,rn,inf);
end
re=prefactor*reshape(re,size(p));
end