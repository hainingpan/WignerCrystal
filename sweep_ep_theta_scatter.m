%Sweep for Wigner Crystal as a function epsilon and theta as scattered
function sweep_ep_theta_scatter(nu,epsilonlist,thetalist)

Ntheta=length(thetalist);
Nep=length(epsilonlist);

load(sprintf('phase1,2_theta(%d,%d,%d).mat',thetalist(1),thetalist(end),Ntheta));

n=27;
final=zeros(Ntheta,Nep);
gap=zeros(Ntheta,Nep);
innergap=zeros(Ntheta,Nep);
finali=zeros(Ntheta,Nep);
reg1=@(x) 76-12*x;
reg2=@(x) 52-9*x;
reg3=@(x) 598/3 - (110*x)/3;
reg4=@(x) 249 - 45*x;
reg5=@(x) 650/3 - (100*x)/3;
reg6=@(x) 305 - 50*x;
for thetai=1:Ntheta
    theta=thetalist(thetai);
    seg1=0:0.1:8;
    seg2=reg1(theta):0.1:reg2(theta);
    seg3=reg3(theta):0.1:reg4(theta);
    seg4=reg5(theta):0.1:reg6(theta);    
    seg=unique(sort([seg1,seg2,seg3,seg4]));
    seg=seg(seg<=60);
    epsilonlist{thetai}=seg;
    thetalist{thetai}=ones(1,length(seg));
end  
    
parfor thetai=1:Ntheta
    parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',thetalist(thetai),'d',10e-9*5.076e6,'nu',nu);
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
    Nepi=length(epsilon(thetai));
    for epi=1:Nepi        
        [final(thetai,epi),spin(:,:,thetai,epi),gap(thetai,epi),innergap(thetai,epi),finali(thetai,epi)]=sweepepsilon(epsilonlist{thetai}(epi),kxlist,kylist,parameters);
    end
end

save(sprintf('phase%d,%d.mat',nu(1),nu(2)),'nu','final','spin','epsilonlist','gap','innergap','thetalist','finali');
end

