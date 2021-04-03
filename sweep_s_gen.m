function sweep_s_gen(nu,theta,slist,hole,perturb,perturbnear,Ush)
% slist=30:10:300;
V1={};
V2={};
energylist={};   
parameters={};
t_bond={};
tlist={};
U_bond={};
Ulist={};
for sindex=1:length(slist)
    s=slist(sindex);
    parameters{sindex}=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'s',s*1e-9*5.076e6,'nu',nu,'Vz',0,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
    if perturb==1
        n=cm(parameters{sindex}.nu(2)/gcd(parameters{sindex}.nu(1),parameters{sindex}.nu(2)),parameters{sindex});
    else
        n=27*(length(parameters{sindex}.Q)<8)+15*(length(parameters{sindex}.Q)>=8)*(length(parameters{sindex}.Q)<16)+9*(length(parameters{sindex}.Q)>=16);
    end
    tshell=3;
    Ushell=length(generate_neighbor(Ush))-1;


    [t,neighborlist]=t_calc_func(tshell,parameters{sindex});
    U=U_calc_func_2(Ushell,parameters{sindex});

    kxlist=zeros(1,n^2);
    kylist=zeros(1,n^2);
    counter=1;
    for xindex=1:n
        for yindex=1:n
            ux=(2*xindex-n-1)/(2*n);
            uy=(2*yindex-n-1)/(2*n);
            klist=ux*parameters{sindex}.bm1+uy*parameters{sindex}.bm2;
            kxlist(counter)=klist(1);
            kylist(counter)=klist(2);
            counter=counter+1;
        end
    end
    kxlist=kxlist';
    kylist=kylist';

    t_bond{sindex}=[neighborlist{1:tshell+1}];
    U_bond{sindex}=[neighborlist{1:Ushell+1}];
    hp=parameters{sindex}.hole;
    tlist{sindex}=-hp*[t{1:tshell+1}];
    Ulist{sindex}=real([U{1:Ushell+1}]);

    parameters{sindex}.N=length(kxlist);
    kxbasis=cell(1,length(parameters{sindex}.Q));
    kybasis=cell(1,length(parameters{sindex}.Q));
    for i=1:length(parameters{sindex}.Q)
        kxbasis{i}=kxlist+parameters{sindex}.Q{i}(1);
        kybasis{i}=kylist+parameters{sindex}.Q{i}(2);
    end

end

parfor sindex=1:length(slist)    
    energylist{sindex}=real(tb(t_bond{sindex},tlist{sindex},[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters{sindex}));

    Qx=cellfun(@(x)x(1),parameters{sindex}.Q);
    Qy=cellfun(@(x)x(2),parameters{sindex}.Q);
    [q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
    [q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
    V1{sindex}=V(U_bond{sindex},Ulist{sindex},q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters{sindex}); %V1_{q_alpha,q_delta}
    [k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
    [k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
    V2{sindex}=V(U_bond{sindex},Ulist{sindex},k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters{sindex}); %V2_{k_alpha,k_beta,q_alpha,q_delta}

end
save(sprintf('phase%d,%d_theta%0.2f_h%d_d%d,%d_U%d.mat',nu(1),nu(2),theta,hole,slist(1),slist(end),Ush),'V1','V2','energylist','slist');
end
