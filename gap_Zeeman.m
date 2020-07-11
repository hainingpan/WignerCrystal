function gap_Zeeman(epsilon,theta)

NVz=51;
Vzlist=linspace(0,50,NVz);
NB=41;
Ezlist=linspace(0,9e-3,NB);

gapmap=zeros(NVz,NB);
szmap=zeros(NVz,NB,3);
isconverge=zeros(NVz,NB);

parfor Vzindex=1:NVz    
parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'Vz',Vzlist(Vzindex),'Ez',0);
tshell=3;
Ushell=0;
[t,neighborlist]=t_calc_func(tshell,parameters);
U=U_calc_func_2(Ushell,parameters);

n=15;
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

t_bond=[neighborlist{1:tshell+1}];
U_bond=[neighborlist{1:Ushell+1}];
hp=1;
tlist=hp*[t{1:tshell+1}];
Ulist=real([U{1:Ushell+1}]);
N=length(kxlist);
kxbasis=cell(1,length(parameters.Q));
kybasis=cell(1,length(parameters.Q));
for i=1:length(parameters.Q)
    kxbasis{i}=kxlist+parameters.Q{i}(1);
    kybasis{i}=kylist+parameters.Q{i}(2);
end
energylist=real(tb(t_bond,tlist,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));
Qx=cellfun(@(x)x(1),parameters.Q);
Qy=cellfun(@(x)x(2),parameters.Q);
[q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
V1=V(U_bond,Ulist,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}
[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
V2=V(U_bond,Ulist,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}
    for Bindex=1:NB
        parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'Vz',Vzlist(Vzindex),'Ez',Ezlist(Bindex));
        parameters.V1=V1;
        parameters.V2=V2;
        parameters.N=N;
        parameters.energylist=energylist;
        [energyall,wfall]=energyMF_init_2(parameters);
        en=zeros(1000,1);
        for i=1:1000
        [en(i),ave,V2deltaave]=totalenergy_2(energyall,wfall,parameters);
        if length(en)>1    
            if abs(en(end)-en(end-1))<1e-15
                break
            end
        end
        [energyall,wfall]=energyMF_2(ave,V2deltaave,parameters);
        end
        [spin,gap,~]=spintexture(energyall,wfall,parameters);
        gapmap(Vzindex,Bindex)=1000*gap;
        isconverge(Vzindex,Bindex)=i;
        szmap(Vzindex,Bindex)=spin(1,4);
    end
end
save(strcat('gap_ep',num2str(epsilon),'theta',num2str(theta),'.mat'),'gapmap','szmap','Vzlist','Ezlist','isconverge');
end