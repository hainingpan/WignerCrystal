Vzlist=-100:10:100;
tstore={};
% Ustore={};
for i=1:length(Vzlist)
    [t]=func(Vzlist(i));
    tstore{i}=t;
%     Ustore{i}=U;
end



function [t]=func(Vz)
parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',[1,2],'d',60e-9*5.076e6,'Vz',Vz,'Ez',0);
tshell=3;
Ushell=length(generate_neighbor(100));
[t,neighborlist]=t_calc_func(tshell,parameters);
% U=U_calc_func_2(Ushell,parameters);
end