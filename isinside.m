function inside=isinside(n,param)
%Check oneNote schematic
%1 if inside 
W=param.W;
L=param.L;
outside=(n(:,2)>=n(:,1)+1)...
    +(n(:,2)<=-n(:,1)-1)...
    +(n(:,2)>=-n(:,1)+2*L)...
    +(n(:,2)<=n(:,1)-2*W);
inside=((outside)==0);
end