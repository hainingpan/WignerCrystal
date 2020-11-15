
% density=[0.6801,0.9740,0.9740,0.9740,0.9928,0.1101,0.1101,0.1101,0.0114,0.0405,0.0114,0.0114];  % nu=[4,24], F, ep=1
% density=[0.6052,0.2281,0.6052,0.2281,0.6052,0.2281,0.2281,0.6052,0.6052,0.6052,0.2281,0.2281]; % nu=[4,24], F, ep=15
% density=[0.2142,0.5032,0.5087,0.5072,0.2139,0.4741,0.4713,0.4756,0.4745,0.2142,0.4702,0.4730]; % nu=[4,24], F, ep=20
% density=[0.4214,0.4120,0.4120,0.4213,0.4214,0.4213,0.4120,0.4120,0.4120,0.4214,0.4120,0.4213]; % nu=[4,24], F, ep=25

density([1,2,4,6,11,3,5,7,8,9,10,12])=density; % permute the index of sites due to the different convention in mainTri2 and mainTMD_2
figure;
hold on;
for i=1:size(param.Aneighbor2,1)
    shift=param.Aneighbor2(i,:)*[param.a1;param.a2];
    scatter(param.r(:,1)+shift(1),param.r(:,2)+shift(2),'k.');
    scatter(param.r(:,1)+shift(1),param.r(:,2)+shift(2),density*40,'r');

%     scatter(param.r(y==1,1)+shift(1),param.r(y==1,2)+shift(2),'r')
end
daspect([1,1,1]);
axis tight;
