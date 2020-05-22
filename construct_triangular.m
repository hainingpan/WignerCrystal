function param=construct_triangular(param)
W=param.W;
L=param.L;
n=zeros((W+1)*(L+1)+W*L,2);
x=0:W+L;
y=-W:L;
[xgrid,ygrid]=meshgrid(x,y);
xline=xgrid(:);
yline=ygrid(:);
filter=isinside([xline,yline],param);
n=[xline(filter),yline(filter)];
param.n=n;
a=[param.a1;param.a2];
param.r=n*a;

nn={};
nn_lin={};
neighbor={[1,-1],[0,-1],[-1,0],[-1,1],[0,1],[1,0]};
for i=1:6
    nn{i}=shift(n+neighbor{i},param);
    nn_lin{i}=arrayfun(@(x) n2i(nn{i}(x,:),n),1:length(n));
end
param.nn_lin=nn_lin;
end