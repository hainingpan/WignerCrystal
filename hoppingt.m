function re=hoppingt(an,parameters)
n=10;
state=1;
xrange=-n:n;
yrange=-n:n;
bM1=parameters.bM1;
bM2=parameters.bM2;

Nkx=length(xrange);
Nky=length(yrange);

%shift to diamond
% a1=-bM1/(2*n);
% a2=(bM1+bM2)/(2*n);
%shift to rectangular
a1=bM2/(2*n);
a2=(-2*bM1-bM2)/2/(2*n);
enmap=zeros(Nkx,Nky);
kxmap=zeros(Nkx,Nky);
kymap=zeros(Nkx,Nky);


omega=abs(cross([bM1,0],[bM2,0]));
omega=omega(3);
parfor xindex=1:Nkx
    kx=xrange(xindex);
    for yindex=1:Nky
        ky=yrange(yindex);
        k=kx*a1+ky*a2;
        [val,~]=energyTMD(k(1),k(2),parameters);
        enmap(xindex,yindex)=val(state);
        kxmap(xindex,yindex)=k(1);
        kymap(xindex,yindex)=k(2);
    end
end

for i=1:length(an)
    neigh=an{i}(1)*parameters.aM1+an{i}(2)*parameters.aM2;

    %use summation
%     re(i)=sum(exp(1i*(kxmap*neigh(1)+kymap*neigh(2))).*enmap,'all')/(Nkx*Nky);

    %use trapzodial integral
    F=exp(1i*(kxmap*neigh(1)+kymap*neigh(2))).*enmap;
    re(i)=trapz(kymap(1,:),trapz(kxmap(:,1),F,2))/omega;
end
end
