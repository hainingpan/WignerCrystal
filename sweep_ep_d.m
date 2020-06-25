%Sweep for Wigner Crystal as a function epsilon and d
function sweep_ep_d(nu)
dlist=linspace(1,10,10);
Nd=length(dlist);
epsilonlist=linspace(1,60,60);
% nu=[1,3];
Nep=length(epsilonlist);
final=zeros(Nd,Nep);
gap=zeros(Nd,Nep);
innergap=zeros(Nd,Nep);

n=15;
kxlist=zeros(1,n^2);
kylist=zeros(1,n^2);
parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',nu);
tshell=3;
Ushell=110;
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

parfor di=1:Nd
param=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'d',dlist(di),'nu',nu);
[t,neighborlist]=t_calc_func(tshell,param);
U=U_calc_func_2(Ushell,param);
sweepfunc=@(x,y) sweepepsilon(tshell,x,y,neighborlist,t,U,kxlist,kylist,param);
    for epi=1:Nep
        [final(di,epi),spin(:,:,di,epi),gap(di,epi),innergap(di,epi)]=sweepfunc(Ushell,epsilonlist(epi));
    end
end
nu=parameters.nu;
save(sprintf('nu%d,%d_U%d.mat',nu(1),nu(2),Ushell),'nu','final','spin','epsilonlist','gap','innergap','dlist');
end

