function nn=shift(n,param)
W=param.W;
L=param.L;
%Check oneNote schematic
shiftvector=(n(:,2)>=n(:,1)+1).*[W,-W]...
    +(n(:,2)<=-n(:,1)-1).*[L,L]...
    +(n(:,2)>=-n(:,1)+2*L).*[-L,-L]...
    +(n(:,2)<=n(:,1)-2*W).*[-W,W];
nn=n+shiftvector;
end