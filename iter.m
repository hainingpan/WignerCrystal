parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',[6,7],'d',10,'Vz',0);
% parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',[1,1],'d',10,'Vz',0);
tshell=3;
Ushell=110;
% [t,neighborlist]=t_calc_func(tshell,parameters);
% U=U_calc_func_2(Ushell,parameters);

n=7;
counter=1;
clear kxlist kylist
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
tlist=-hp*[t{1:tshell+1}];
epsilon=10;
Ulist=real([U{1:Ushell+1}])/epsilon;

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
% parameters.V1=V1;
% parameters.V1=zeros(length(Qx),length(Qy));
[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
parameters.V2=V(U_bond,Ulist,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}
% parameters.V2=V2;
% parameters.V2=zeros(length(kxlist),length(kylist),length(Qx),length(Qy));


clear spinsav en gapsav
[energyall,wfall]=energyMF_init_2(parameters);
for i=1:1000
[spin,gap]=spintexture(energyall,wfall,parameters);
en(i)=totalenergy_2(energyall,wfall,parameters);
fprintf("%d: gap:%0.8f meV E:%f meV\n",i,1000*gap,1000*en(end));
disp([spin,angle(spin(:,2)+spin(:,3)*1i)*180/pi,angle(spin(:,4)+sqrt(spin(:,3).^2+spin(:,2).^2)*1i)*180/pi])
plot(en);
ylabel('energy (eV)');
drawnow;
spinsav(:,:,i)=spin;
gapsav(i)=gap;
[energyall,wfall]=energyMF_2(energyall,wfall,parameters);
if length(en)>1    
    if abs(en(end)-en(end-1))<1e-12
        break
    end
end
end
final=en(end);
% save(sprintf('nu%d,%d_t%d_U%d_hp%d_ep%d.mat',parameters.nu(1),parameters.nu(2),tshell,Ushell,hp,epsilon),'en','spinsav','gapsav')
