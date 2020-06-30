function [Elist,dos]=DOS(energyall,parameters)
energyall_sort=sort(energyall(:));
mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));

eta=1e-4;
Emin=min(energyall,[],'all');
Emax=max(energyall,[],'all');
Elist=linspace(Emin,Emax,1000);
energyproj=energyall(:);
deltaf=1/pi*eta./((Elist-energyproj).^2+eta^2); %axis 1: state; axis 2: Elist

dos=sum(deltaf,1);

dos=dos/size(energyall,1)/(sqrt(3)/2*(parameters.aM/5.076e-3)^2);  % of eV^-1 nm^-2
% Elist=Elist-mu;
% save(sprintf('dos_nu%d,%d.mat',parameters.nu(1),parameters.nu(2)),'Elist','dos');
% figure;plot(1000*Elist,dos);axis tight;
% savefig(sprintf('dos_nu%d,%d.fig',parameters.nu(1),parameters.nu(2)));
nlist=arrayfun(@(x) sum(x>energyall_sort),Elist)/length(energyall_sort);
figure;plot(nlist,dos);axis tight;
end

