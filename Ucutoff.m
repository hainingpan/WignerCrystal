% Ur2=13.6./(rlist*3.28e-10/parameters.theta/5.29177210903e-11)*2...
% -13.6./(sqrt((rlist*3.28e-10/parameters.theta).^2+(1*3.28e-10/parameters.theta).^2)/5.29177210903e-11)*2;
Ur=cellfun(@(x)real(x(1)),U);
neighbor=arrayfun(@(i) norm(neighborlist{i}{1}*[0,-1;sqrt(3)/2,-1/2]),1:length(U));
alpha=0.00729735;
d=parameters.d;
rlist=0.5:0.1:20;
Ur2=alpha./(rlist*parameters.aM)-alpha./(sqrt((rlist*parameters.aM).^2+parameters.d.^2));

figure;
plot(rlist,Ur2);
hold on;
plot(neighbor(1:length(Ur)),Ur(1:end),'.');
xlim([0,20])