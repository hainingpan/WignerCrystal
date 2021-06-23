function sweep_vz_gen(nu,theta,vzlist,hole,perturb,perturbnear,Ush)
V1={};
V2={};
energylist={};   
parameters={};
t_bond={};
tlist={};
U_bond={};
Ulist={};
for vzindex=1:length(vzlist)
    vz=vzlist(vzindex);
    parameters{vzindex}=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'d',inf,'nu',nu,'Vz',0,'hole',hole,'perturb',perturb,'perturbnear',perturbnear,'vz',vz);
    if perturb==1
        n=cm(parameters{vzindex}.nu(2)/gcd(parameters{vzindex}.nu(1),parameters{vzindex}.nu(2)),parameters{vzindex});
    else
        n=27*(length(parameters{vzindex}.Q)<8)+15*(length(parameters{vzindex}.Q)>=8)*(length(parameters{vzindex}.Q)<16)+9*(length(parameters{vzindex}.Q)>=16);
    end
    tshell=3;
    Ushell=length(generate_neighbor(Ush))-1;


    [t,neighborlist]=t_calc_func(tshell,parameters{vzindex});
    U=U_calc_func_2(Ushell,parameters{vzindex});

    kxlist=zeros(1,n^2);
    kylist=zeros(1,n^2);
    counter=1;
    for xindex=1:n
        for yindex=1:n
            ux=(2*xindex-n-1)/(2*n);
            uy=(2*yindex-n-1)/(2*n);
            klist=ux*parameters{vzindex}.bm1+uy*parameters{vzindex}.bm2;
            kxlist(counter)=klist(1);
            kylist(counter)=klist(2);
            counter=counter+1;
        end
    end
    kxlist=kxlist';
    kylist=kylist';

    t_bond{vzindex}=[neighborlist{1:tshell+1}];
    U_bond{vzindex}=[neighborlist{1:Ushell+1}];
    hp=parameters{vzindex}.hole;
    tlist{vzindex}=-hp*[t{1:tshell+1}];
    Ulist{vzindex}=real([U{1:Ushell+1}]);

    parameters{vzindex}.N=length(kxlist);
    kxbasis=cell(1,length(parameters{vzindex}.Q));
    kybasis=cell(1,length(parameters{vzindex}.Q));
    for i=1:length(parameters{vzindex}.Q)
        kxbasis{i}=kxlist+parameters{vzindex}.Q{i}(1);
        kybasis{i}=kylist+parameters{vzindex}.Q{i}(2);
    end

end

parfor vzindex=1:length(vzlist)    
    energylist{vzindex}=real(tb(t_bond{vzindex},tlist{vzindex},[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters{vzindex}));

    Qx=cellfun(@(x)x(1),parameters{vzindex}.Q);
    Qy=cellfun(@(x)x(2),parameters{vzindex}.Q);
    [q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
    [q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
    V1{vzindex}=V(U_bond{vzindex},Ulist{vzindex},q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters{vzindex}); %V1_{q_alpha,q_delta}
    [k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
    [k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
    V2{vzindex}=V(U_bond{vzindex},Ulist{vzindex},k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters{vzindex}); %V2_{k_alpha,k_beta,q_alpha,q_delta}

end
save(sprintf('phase%d,%d_theta%0.2f_h%d_vz%d,%d_U%d.mat',nu(1),nu(2),theta,hole,vzlist(1),vzlist(end),Ush),'V1','V2','energylist','vzlist','t','U','-v7.3');
end
