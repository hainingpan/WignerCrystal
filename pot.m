function sum=pot(x,y,rank,param)
sum=0;
X=param.L*(param.a1+param.a2);
Y=param.W*(param.a1-param.a2);
for i=-rank:rank
    for j=-rank:rank
        y2=y+i*X+j*Y;
        r=dist(x,y2,param);
        en=1/r-1/sqrt(r^2+param.d^2);
        sum=sum+en;
    end
end

end