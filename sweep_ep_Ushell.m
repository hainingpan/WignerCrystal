%Sweep for Wigner Crystal as a function epsilon and theta
function sweep_ep_Ushell(nu,epsilonlist,U_sh_list,hole,perturb,perturbnear,filename,theta)

NU_sh_list=length(U_sh_list);
Nep=length(epsilonlist);
energylist={};
V1={};
V2={};
load(filename);
param=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'d',60e-9*5.076e6,'nu',nu,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
if perturb==1
%     n=cm(abs(1/(1-nu(1)/nu(2))));
    n=cm(param.nu(2)/gcd(param.nu(1),param.nu(2)),param);
else
    n=27*(length(param.Q)<8)+15*(length(param.Q)>=8)*(length(param.Q)<16)+9*(length(param.Q)>=16);
end
final=zeros(NU_sh_list,Nep);
gap=zeros(NU_sh_list,Nep);
innergap=zeros(NU_sh_list,Nep);
finali=zeros(NU_sh_list,Nep);
ch=zeros(NU_sh_list,Nep);
parfor U_shi=1:NU_sh_list
    parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'d',60e-9*5.076e6,'nu',nu,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
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

    parameters.N=length(kxlist);
    kxbasis=cell(1,length(parameters.Q));
    kybasis=cell(1,length(parameters.Q));
    for i=1:length(parameters.Q)
        kxbasis{i}=kxlist+parameters.Q{i}(1);
        kybasis{i}=kylist+parameters.Q{i}(2);
    end
    parameters.energylist=energylist{U_shi};

    parameters.V1=V1{U_shi};
    parameters.V2=V2{U_shi};
    for epi=1:Nep
        [final(U_shi,epi),spin(:,:,U_shi,epi),gap(U_shi,epi),innergap(U_shi,epi),finali(U_shi,epi),ch(U_shi,epi)]=sweepepsilon(epsilonlist(epi),kxlist,kylist,parameters);
    end
end
save(sprintf('phase%d,%d_h%d_theta%0.2f.mat',nu(1),nu(2),hole,theta),'nu','final','spin','epsilonlist','gap','innergap','U_sh_list','finali','ch');

end

