%Sweep for Wigner Crystal as a function epsilon and theta
function tU_ep_theta_gen(thetalist)

Ntheta=length(thetalist);

tshell=3;
Ushell=0;
Ushell=length(generate_neighbor(100));

t={};
U={};
for thetai=1:Ntheta
    disp(thetai);
    parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',thetalist(thetai),'d',60e-9*5.076e6,'nu',[1,2]);
    [t{thetai},neighborlist]=t_calc_func(tshell,parameters);
    U{thetai}=U_calc_func_2(Ushell,parameters);    
end
save(sprintf('tU_theta(%0.2f,%0.2f,%d).mat',thetalist(1),thetalist(end),Ntheta),'t','U');
end

