function sweep_d_gen(nu,theta,dlist,hole,perturb,perturbnear)
% dlist=30:10:300;
V1={};
V2={};
energylist={};   
parameters={};
t_bond={};
tlist={};
U_bond={};
Ulist={};
for dindex=1:length(dlist)
    d=dlist(dindex);
    parameters{dindex}=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'d',d*1e-9*5.076e6,'nu',nu,'Vz',0,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
    if perturb==1
        n=cm(parameters{dindex}.nu(2)/gcd(parameters{dindex}.nu(1),parameters{dindex}.nu(2)),parameters{dindex});
    else
        n=27*(length(parameters{dindex}.Q)<8)+15*(length(parameters{dindex}.Q)>=8)*(length(parameters{dindex}.Q)<16)+9*(length(parameters{dindex}.Q)>=16);
    end
    tshell=3;
    Ushell=length(generate_neighbor(100))-1;


    [t,neighborlist]=t_calc_func(tshell,parameters{dindex});
    U=U_calc_func_2(Ushell,parameters{dindex});

    kxlist=zeros(1,n^2);
    kylist=zeros(1,n^2);
    counter=1;
    for xindex=1:n
        for yindex=1:n
            ux=(2*xindex-n-1)/(2*n);
            uy=(2*yindex-n-1)/(2*n);
            klist=ux*parameters{dindex}.bm1+uy*parameters{dindex}.bm2;
            kxlist(counter)=klist(1);
            kylist(counter)=klist(2);
            counter=counter+1;
        end
    end
    kxlist=kxlist';
    kylist=kylist';

    t_bond{dindex}=[neighborlist{1:tshell+1}];
    U_bond{dindex}=[neighborlist{1:Ushell+1}];
    hp=parameters{dindex}.hole;
    tlist{dindex}=-hp*[t{1:tshell+1}];
    Ulist{dindex}=real([U{1:Ushell+1}]);

    parameters{dindex}.N=length(kxlist);
    kxbasis=cell(1,length(parameters{dindex}.Q));
    kybasis=cell(1,length(parameters{dindex}.Q));
    for i=1:length(parameters{dindex}.Q)
        kxbasis{i}=kxlist+parameters{dindex}.Q{i}(1);
        kybasis{i}=kylist+parameters{dindex}.Q{i}(2);
    end

end

parfor dindex=1:length(dlist)    
    energylist{dindex}=real(tb(t_bond{dindex},tlist{dindex},[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters{dindex}));

    Qx=cellfun(@(x)x(1),parameters{dindex}.Q);
    Qy=cellfun(@(x)x(2),parameters{dindex}.Q);
    [q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
    [q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
    V1{dindex}=V(U_bond{dindex},Ulist{dindex},q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters{dindex}); %V1_{q_alpha,q_delta}
    [k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
    [k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
    V2{dindex}=V(U_bond{dindex},Ulist{dindex},k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters{dindex}); %V2_{k_alpha,k_beta,q_alpha,q_delta}

end
save(sprintf('phase%d,%d_theta%0.2f_h%d_d%d,%d.mat',nu(1),nu(2),theta,hole,dlist(1),dlist(-1)),'V1','V2','energylist','dlist');
end
