clear
thetalist=linspace(3,5,10);
% masslist=linspace(0.1,1,10);
loclist=thetalist*0;
for index=1:length(thetalist)
    disp(index)
parameters=mainTMD_2('m',0.65,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',thetalist(index),'nu',[1,1]*1,...
    'd',60e-9*5.076e6,'Vz',0,'Ez',0,'hole',1,'perturb',0,'perturbnear',[1,3]);
loclist(index)=wannier_localization_length(parameters);
end
storelist=[thetalist,loclist];

save('theta_mass0.65.mat')
