%Sweep for Wigner Crystal as a function epsilon and theta
function sweep_ep_theta(nu,epsilonlist,thetalist,filename)

Ntheta=length(thetalist);
Nep=length(epsilonlist);

% filename='phase4,12_theta(3.00,5.00,51).mat';
load(filename);
% load(sprintf('phase1,2_theta(%.2f,%.2f,%d)_d60.mat',thetalist(1),thetalist(end),Ntheta));

n=27*(length(param.Q)<=8)+15*(length(param.Q)>8;
final=zeros(Ntheta,Nep);
gap=zeros(Ntheta,Nep);
innergap=zeros(Ntheta,Nep);
finali=zeros(Ntheta,Nep);
ch=zeros(Ntheta,Nep);
parfor thetai=1:Ntheta
    parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',thetalist(thetai),'d',60e-9*5.076e6,'nu',nu);
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
    parameters.energylist=energylist{thetai};

    parameters.V1=V1{thetai};
    parameters.V2=V2{thetai};
    for epi=1:Nep
        [final(thetai,epi),spin(:,:,thetai,epi),gap(thetai,epi),innergap(thetai,epi),finali(thetai,epi),ch(thetai,epi)]=sweepepsilon(epsilonlist(epi),kxlist,kylist,parameters);
    end
end
save(sprintf('phase%d,%d.mat',nu(1),nu(2)),'nu','final','spin','epsilonlist','gap','innergap','thetalist','finali','ch');
end

