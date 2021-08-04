function sweep_vz_T(nu,vzlist,Tlist,theta,hole,perturb,perturbnear,filename,epsilon,Ush)

    NT=length(Tlist);
    Nvz=length(vzlist);
    energylist={};
    % thetalist=linspace(3,5,51);
    % assert(sum(thetalist==theta)==1,'Theta=%.2f not found',theta)
    % thetai=find(thetalist==theta);
    V1={};
    V2={};
    load(filename);
    param=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'T',0,'d',inf*5.076e6,'nu',nu,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
    if perturb==1
        n=cm(param.nu(2)/gcd(param.nu(1),param.nu(2)),param);
    else
        n=27*(length(param.Q)<8)+15*(length(param.Q)>=8)*(length(param.Q)<16)+9*(length(param.Q)>=16);
    end
    final=zeros(NT,Nvz);
    gap=zeros(NT,Nvz);
    innergap=zeros(NT,Nvz);
    finali=zeros(NT,Nvz);
    ch=zeros(NT,Nvz);
    parfor Ti=1:NT
        for vzi=1:Nvz
            parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'T',Tlist(Ti),'d',inf,'vz',vzlist(vzi),'nu',nu,'hole',hole,'perturb',perturb,'perturbnear',perturbnear);
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
            parameters.energylist=energylist{vzi};
    
            parameters.V1=V1{vzi};
            parameters.V2=V2{vzi};
            [final(Ti,vzi),spin(:,:,Ti,vzi),gap(Ti,vzi),innergap(Ti,vzi),finali(Ti,vzi),ch(Ti,vzi)]=sweepd(epsilon,kxlist,kylist,parameters);
        end
    end
    save(sprintf('phase%d,%d_h%d_vz_ep%d_theta%.2f_U%d.mat',nu(1),nu(2),hole,epsilon,theta,Ush),'nu','final','spin','vzlist','gap','innergap','Tlist','finali','ch');
    
    end
    
    