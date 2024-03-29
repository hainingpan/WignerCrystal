function sweep_Ushell_gen(nu,theta,hole,perturb,perturbnear,U_sh)

parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'d',60e-9*5.076e6,'nu',nu,'Vz',0,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
if perturb==1
%     n=cm(abs(1/(1-nu(1)/nu(2))));
    n=cm(parameters.nu(2)/gcd(parameters.nu(1),parameters.nu(2)),parameters);
else
    n=27*(length(parameters.Q)<8)+15*(length(parameters.Q)>=8)*(length(parameters.Q)<16)+9*(length(parameters.Q)>=16);
end
tshell=3;
Ushell=length(generate_neighbor(U_sh))-1;
V1={};
V2={};
energylist={};

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
hp=parameters.hole;
tlist=-hp*[t{1:tshell+1}];
Ulist=real([U{1:Ushell+1}]);

parameters.N=length(kxlist);
kxbasis=cell(1,length(parameters.Q));
kybasis=cell(1,length(parameters.Q));
for i=1:length(parameters.Q)
    kxbasis{i}=kxlist+parameters.Q{i}(1);
    kybasis{i}=kylist+parameters.Q{i}(2);
end
energylist=real(tb(t_bond,tlist,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));

Qx=cellfun(@(x)x(1),parameters.Q);
Qy=cellfun(@(x)x(2),parameters.Q);
[q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
V1=V(U_bond,Ulist,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}
[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
V2=V(U_bond,Ulist,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}

save(sprintf('phase%d,%d_theta%0.2f_h%d_U%d.mat',nu(1),nu(2),theta,hole,U_sh),'V1','V2','energylist','t','U');
end
