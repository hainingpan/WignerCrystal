[kxq,kyq]=meshgrid(linspace(min(kxlist),max(kxlist),100),linspace(min(kylist),max(kylist),100));
energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
figure;
hold on;
surf(2*kxq/norm(parameters.bm1),2*kyq/norm(parameters.bm2),1000*mu*ones(100),'edgecolor','none','FaceAlpha',0.2);
% for i=1:size(energyall,2)
for i=2:3
    vq=griddata(kxlist,kylist,1000*energyall(:,i),kxq,kyq);
    mesh(kxq/norm(parameters.bm1),kyq/norm(parameters.bm1),vq);
end
xlabel('k_x/|b_M|');
ylabel('k_y/|b_M|');
zlabel('E (meV)');
title(sprintf("%d bands\ngap: %0.8f meV",size(energyall,2),1000*gap));
xlim([min(kxq(:)),max(kxq(:))]/norm(parameters.bm1));
ylim([min(kyq(:)),max(kyq(:))]/norm(parameters.bm1));

view([-1,-1,.1]);
