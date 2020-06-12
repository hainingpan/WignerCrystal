function sum=pot2(x,y,param)
sum=0;
for i=1:size(param.Aneighbor2,1)
    y2=y+param.Aneighbor2(i,:);
    r=norm((y2-x)*[param.a1;param.a2]);
    en=(1/r-1/sqrt(r^2+param.d^2))*(r<=param.neighbor(2)*1.01);
    sum=sum+en;
end

end