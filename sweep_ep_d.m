%Sweep for Wigner Crystal as a function epsilon and d
function sweep_ep_d(nu,epsilonlist,dlist)
% tic;
% dlist=linspace(1,10,20);
Nd=length(dlist);
% epsilonlist=linspace(1,80,160);
% nu=[1,3];
Nep=length(epsilonlist);
final=zeros(Nd,Nep);
gap=zeros(Nd,Nep);
innergap=zeros(Nd,Nep);

param=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',nu);
n=21*(length(param.Q)<=8)+15*(length(param.Q)>8 & length(param.Q)<=21)+9*(length(param.Q)>21);
kxlist=zeros(1,n^2);
kylist=zeros(1,n^2);
tshell=3;
Ushell=373;
counter=1;
for xindex=1:n
    for yindex=1:n
        ux=(2*xindex-n-1)/(2*n);
        uy=(2*yindex-n-1)/(2*n);
        klist=ux*param.bm1+uy*param.bm2;
        kxlist(counter)=klist(1);
        kylist(counter)=klist(2);
        counter=counter+1;
    end
end
kxlist=kxlist';
kylist=kylist';

final=zeros(Nd,Nep);
gap=zeros(Nd,Nep);
innergap=zeros(Nd,Nep);
finali=zeros(Nd,Nep);
parfor di=1:Nd
    parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'d',dlist(di),'nu',nu);
    [t,neighborlist]=t_calc_func(tshell,parameters);
    U=U_calc_func_2(Ushell,parameters);

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
    parameters.energylist=real(tb(t_bond,tlist,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));

    Qx=cellfun(@(x)x(1),parameters.Q);
    Qy=cellfun(@(x)x(2),parameters.Q);
    [q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
    [q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
    parameters.V1=V(U_bond,Ulist,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}
    [k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
    [k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
    parameters.V2=V(U_bond,Ulist,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}
    for epi=1:Nep
        [final(di,epi),spin(:,:,di,epi),gap(di,epi),innergap(di,epi),finali(di,epi)]=sweepepsilon(epsilonlist(epi),parameters);
    end
end
save(sprintf('nu%d,%d_U%d.mat',param.nu(1),param.nu(2),Ushell),'nu','final','spin','epsilonlist','gap','innergap','dlist','finali');
end

