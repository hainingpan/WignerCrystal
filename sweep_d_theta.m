function sweep_d_theta(nu,dlist,thetalist,hole,perturb,perturbnear,filename,epsilon,Ush)

Ntheta=length(thetalist);
Nd=length(dlist);
energylist={};
V1={};
V2={};
% filename='phase4,12_theta(3.00,5.00,51).mat';
load(filename);
% load(sprintf('phase1,2_theta(%.2f,%.2f,%d)_d60.mat',thetalist(1),thetalist(end),Ntheta));
param=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',3,'d',60e-9*5.076e6,'nu',nu,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
if perturb==1
%     n=cm(abs(1/(1-nu(1)/nu(2))));
    n=cm(param.nu(2)/gcd(param.nu(1),param.nu(2)),param);
else
    n=27*(length(param.Q)<8)+15*(length(param.Q)>=8)*(length(param.Q)<16)+9*(length(param.Q)>=16);
end
final=zeros(Ntheta,Nd);
gap=zeros(Ntheta,Nd);
innergap=zeros(Ntheta,Nd);
finali=zeros(Ntheta,Nd);
ch=zeros(Ntheta,Nd);
parfor thetai=1:Ntheta
    for di=1:Nd
        parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',thetalist(thetai),'d',dlist(di)*1e-9*5.076e6,'nu',nu,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
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
        parameters.energylist=energylist{thetai}{di};

        parameters.V1=V1{thetai}{di};
        parameters.V2=V2{thetai}{di};
        [final(thetai,di),spin(:,:,thetai,di),gap(thetai,di),innergap(thetai,di),finali(thetai,di),ch(thetai,di)]=sweepd(epsilon,kxlist,kylist,parameters);
    end
end
save(sprintf('phase%d,%d_h%d_d_ep%d_U%d.mat',nu(1),nu(2),hole,epsilon,Ush),'nu','final','spin','dlist','gap','innergap','thetalist','finali','ch');

end

