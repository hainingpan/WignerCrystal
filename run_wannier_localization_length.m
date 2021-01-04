function run_wannier_localization_length(m)
thetalist=linspace(1,5,20);
% masslist=linspace(0.1,1,10);
loclist=thetalist*0;
rs_inv=thetalist*0;
tlist=thetalist*0;
Ulist=thetalist*0;
for index=1:length(thetalist)
    disp(index)
    parameters=mainTMD_2('m',m,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',thetalist(index),'nu',[2,3]*1,...
        'd',60e-9*5.076e6,'Vz',0,'Ez',0,'hole',1,'perturb',0,'perturbnear',[1,3]);
    loclist(index)=wannier_localization_length(parameters);
    rs_inv(index)=rs(1,parameters);

    [t,neighborlist]=t_calc_func(1,parameters);
    U=U_calc_func_2(0,parameters);
    tlist(index)=real(t{2}(1));
    Ulist(index)=U{1};
end
% storelist=[thetalist,loclist];
save(strcat('theta_mass_',num2str(m),'_nu2,3.mat'));
end
