%Sweep for Wigner Crystal as a function epsilon and theta
function sweep_ep_theta_gen(nu,epsilonlist,thetalist,Vz)

Ntheta=length(thetalist);
Nep=length(epsilonlist);

% n=27;
param=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',3,'d',60e-9*5.076e6,'nu',nu,'Vz',Vz);
n=27*(length(param.Q)<8)+15*(length(param.Q)>=8);
tshell=3;
Ushell=length(generate_neighbor(100));
V1={};
V2={};
energylist={};
for thetai=1:Ntheta
    disp(thetai);
    parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',thetalist(thetai),'d',60e-9*5.076e6,'nu',nu,'Vz',Vz);
    [t,neighborlist]=t_calc_func(tshell,parameters);
    U=U_calc_func_2(Ushell,parameters);
    
    kxlist=zeros(1,n^2);
    kylist=zeros(1,n^2);
    counter=1;
    for xindex=1:n
        for yindex=1:n
            ux=(2*xindex-n-1)/(2*n);
            uy=(2*yindex-n-1)/(2*n);
            klist=ux*parameters.bm1+uy*parameters.bm2;
            kxlist(counter)=klist(1);
            kylist(counter)=klist(2);
            counter=counter+1;
        end
    end
    kxlist=kxlist';
    kylist=kylist';

    t_bond=[neighborlist{1:tshell+1}];
    U_bond=[neighborlist{1:Ushell+1}];
    hp=1;
    tlist=-hp*[t{1:tshell+1}];
    Ulist=real([U{1:Ushell+1}]);

    parameters.N=length(kxlist);
    kxbasis=cell(1,length(parameters.Q));
    kybasis=cell(1,length(parameters.Q));
    for i=1:length(parameters.Q)
        kxbasis{i}=kxlist+parameters.Q{i}(1);
        kybasis{i}=kylist+parameters.Q{i}(2);
    end
    energylist{thetai}=real(tb(t_bond,tlist,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));

    Qx=cellfun(@(x)x(1),parameters.Q);
    Qy=cellfun(@(x)x(2),parameters.Q);
    [q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
    [q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
    V1{thetai}=V(U_bond,Ulist,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}
    [k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
    [k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
    V2{thetai}=V(U_bond,Ulist,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}

end
save(sprintf('phase%d,%d_theta(%0.2f,%0.2f,%d)_Vz(%0.1f).mat',nu(1),nu(2),thetalist(1),thetalist(end),Ntheta,Vz),'V1','V2','energylist','t','U');
end

