Ur=cellfun(@(x)x(1),U);
% Ushell=373;
nlist=cellfun(@(x) norm(x{1}(1)*parameters.aM1+x{1}(2)*parameters.aM2),neighborlist(1:Ushell+1),'UniformOutput',false);
nlist=cell2mat(nlist);
alpha=0.00729735;
Ur2=alpha./(nlist)-alpha./sqrt(nlist.^2+(parameters.d)^2);
figure;plot(nlist/parameters.aM,Ur(1:Ushell+1),nlist/parameters.aM,Ur2);
legend('origin','fit')
