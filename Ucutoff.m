% Ur2=13.6./(rlist*3.28e-10/parameters.theta/5.29177210903e-11)*2...
% -13.6./(sqrt((rlist*3.28e-10/parameters.theta).^2+(1*3.28e-10/parameters.theta).^2)/5.29177210903e-11)*2;
Ur=cellfun(@(x)x(1),U);

alpha=0.00729735;
d=parameters.d;
rlist=1:0.1:9;
Ur2=alpha./(rlist*parameters.aM)-alpha./(sqrt((rlist*parameters.aM).^2+parameters.d.^2));

figure;plot(neighbor(2:length(Ur)),Ur(2:end),rlist,Ur2)