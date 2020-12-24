parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',[1,1]*1,...
    'd',60e-9*5.076e6,'Vz',0,'Ez',0,'hole',1,'perturb',0,'perturbnear',[1,3]);
tshell=3;
% tshell=1;

Ushell=length(generate_neighbor(100));
% % Ushell=0;
[t,neighborlist]=t_calc_func(tshell,parameters);
U=U_calc_func_2(Ushell,parameters);

t=cellfun(@(x) mean(x)*ones(1,length(x)),t,'UniformOutput',false);

epsilon=1;

n=15;
% n=cm(abs(1/(1-parameters.nu(1)/parameters.nu(2))),parameters);
% n=cm(parameters.nu(2)/gcd(parameters.nu(1),parameters.nu(2)),parameters);

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
hp=parameters.hole;
tlist=-hp*[t{1:tshell+1}];  %(-) is for hole- like doping; (+) is for electron-like doping
Ulist=real([U{1:Ushell+1}])/1;

parameters.N=length(kxlist);
kxbasis=cell(1,length(parameters.Q));
kybasis=cell(1,length(parameters.Q));
for i=1:length(parameters.Q)
    kxbasis{i}=kxlist+parameters.Q{i}(1);
    kybasis{i}=kylist+parameters.Q{i}(2);
end
energylist=real(tb(t_bond,tlist,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));
parameters.energylist=energylist;

Qx=cellfun(@(x)x(1),parameters.Q);
Qy=cellfun(@(x)x(2),parameters.Q);
[q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
[q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
V1=V(U_bond,Ulist,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}
parameters.V1=V1/epsilon;
% parameters.V1=zeros(length(Qx),length(Qy));
[k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
[k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
V2=V(U_bond,Ulist,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}
parameters.V2=V2/epsilon;
% parameters.V2=zeros(length(kxlist),length(kylist),length(Qx),length(Qy));

clear spinsav en gapsav chsav
if ismember(parameters.nu,[[4,8];[5,10];[2,2];[4,4];[7,28];[8,32];[9,12];[21,28];[24,32]],'row')
    [ave,V2ave]=average_kagome(parameters.phi,parameters.s1,kxlist,kylist,parameters);
    [energyall,wfall]=energyMF_2(ave,V2ave,parameters);
else
    if parameters.nu==[4,12]
        [ave,V2ave]=average_honeycomb(kxlist,kylist,parameters);
        [energyall,wfall]=energyMF_2(ave,V2ave,parameters);
    else
        if ismember(parameters.nu,[[4,6];[6,6]],'row')
            [ave,V2ave]=average_honeycomb2(kxlist,kylist,parameters);
            [energyall,wfall]=energyMF_2(ave,V2ave,parameters);        
        else
            [energyall,wfall]=energyMF_init_2(parameters);
        end
    end
end

fig1=figure;
% fig2=figure;
for i=1:300
[spin,gap,innergap]=spintexture(energyall,wfall,parameters);
[en(i),ave,V2deltaave]=totalenergy_2(energyall,wfall,parameters);
chsav(i)=chern(wfall,parameters);
fprintf("%d: gap:%0.8f meV E:%f meV innergap: %0.8f\n",i,1000*gap,1000*en(end),1000*innergap);
disp([spin,angle(spin(:,2)+spin(:,3)*1i)*180/pi,angle(spin(:,4)+sqrt(spin(:,3).^2+spin(:,2).^2)*1i)*180/pi])
% figure(fig1);
plot(en);
ylabel('energy (eV)');
% plotband;
drawnow;
spinsav(:,:,i)=spin;
% plot(squeeze(angle(spinsav(1,4,:)+sqrt(spinsav(1,3,:).^2+spinsav(1,2,:).^2)*1i)*180/pi));
gapsav(i)=gap;
if length(en)>1     
    if abs(en(end)-en(end-1))<1e-15
        break
    end
end
[energyall,wfall]=energyMF_2(ave,V2deltaave,parameters);
end
final=en(end);
plotband;
% save(sprintf('nu%d,%d_t%d_U%d_hp%d_ep%d.mat',parameters.nu(1),parameters.nu(2),tshell,Ushell,hp,epsilon),'en','spinsav','gapsav')
