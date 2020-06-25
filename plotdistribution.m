% figure;scatter(param.r(:,1),param.r(:,2))
% y=[1 1 0 0 1 0 0 0 1 1 1 0]; %rectangular
% y=[0,1,1,1,1,1,1,0,0,0,0,0]; %hexagonal
% y=[1 1 0 1 1 0 1 1 1 1 1 0]; %kagome
% y=[1,1,1,1,1,1,1,0,1,0,1,0]; %complementary kagome
% y=[0,1,1,1]; %complementary rectangular
% y=~[0,1,1,1]; %rectangular
% y=[0 1 0 0 0 1 0 0 0]; %honeycomb 1/9
% y=[1 0 1 0 1 0 1 0 0]; % 4/9
% y=[0 0 1 0 1 0 1 0 0]; %3/9
% y=[0 0 0 1 0 0 1]; %2/7
% y=[1 0 0 1 0 0 1]; %3/7
% y=[1 1 1 0 1 1 0]; %5/7
%%
% y=[1 0 0 0 0]; %[1,5]
% y=[0 1 1 0 0]; %[2,5]
% y=[0 1 1 0 1]; %[3,5]
y=~[1 0 0 0 0]; %[4,5]
%% complete
figure;
hold on;
for i=1:size(param.Aneighbor2,1)
    shift=param.Aneighbor2(i,:)*[param.a1;param.a2];
    scatter(param.r(:,1)+shift(1),param.r(:,2)+shift(2),'k.');
    scatter(param.r(y==1,1)+shift(1),param.r(y==1,2)+shift(2),'r')
end

%% unit cell
figure;scatter(param.r(:,1),param.r(:,2))
hold on;scatter(param.r(y==1,1),param.r(y==1,2),'.')