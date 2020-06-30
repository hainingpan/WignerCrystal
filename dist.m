function re=dist(x,y,param)
% re=norm(x-y);
rel=x-y;
L=param.L*sqrt(3);
W=param.W;
re=sqrt((L/2-mod(rel(1)+L/2,L))^2+(W/2-mod(rel(2)+W/2,W))^2);
end