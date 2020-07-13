% parameters=mainTMD_spin('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',3,'nu',[1,2],'d',10);
% tshell=3;
% Ushell=length(generate_neighbor(100));
% [t,neighborlist]=t_calc_func(tshell,parameters);
% U=U_calc_func_2(Ushell,parameters);
% 
% n=15;
% counter=1;
% clear kxlist kylist
% for xindex=1:n
%     for yindex=1:n
%         ux=(2*xindex-n-1)/(2*n);
%         uy=(2*yindex-n-1)/(2*n);
%         klist=ux*parameters.bm1+uy*parameters.bm2;
%         kxlist(counter)=klist(1);
%         kylist(counter)=klist(2);
%         counter=counter+1;
%     end
% end
% kxlist=kxlist';
% kylist=kylist';
% 
% t_bond=[neighborlist{1:tshell+1}];
% U_bond=[neighborlist{1:Ushell+1}];
% hp=1;
% tlist=-hp*[t{1:tshell+1}];
epsilon=30;
% Ulist=real([U{1:Ushell+1}])/1;
% % 
% parameters.N=length(kxlist);
% kxbasis=cell(1,length(parameters.Q));
% kybasis=cell(1,length(parameters.Q));
% for i=1:length(parameters.Q)
%     kxbasis{i}=kxlist+parameters.Q{i}(1);
%     kybasis{i}=kylist+parameters.Q{i}(2);
% end
% parameters.energylist=real(tb(t_bond,tlist,[cell2mat(kxbasis),-cell2mat(kxbasis)],[cell2mat(kybasis),-cell2mat(kybasis)],parameters));
% 
% Qx=cellfun(@(x)x(1),parameters.Q);
% Qy=cellfun(@(x)x(2),parameters.Q);
% [q_alpha_x,q_delta_x]=meshgrid(Qx,Qx);
% [q_alpha_y,q_delta_y]=meshgrid(Qy,Qy);
% V1=V(U_bond,Ulist,q_alpha_x-q_delta_x,q_alpha_y-q_delta_y,parameters); %V1_{q_alpha,q_delta}
parameters.V1=V1/epsilon;
% [k_alpha_x,k_beta_x,q_alpha_x,q_delta_x]=ndgrid(kxlist,kxlist,Qx,Qx);
% [k_alpha_y,k_beta_y,q_alpha_y,q_delta_y]=ndgrid(kylist,kylist,Qy,Qy);
% V2=V(U_bond,Ulist,k_alpha_x-k_beta_x+q_alpha_x-q_delta_x,k_alpha_y-k_beta_y+q_alpha_y-q_delta_y,parameters); %V2_{k_alpha,k_beta,q_alpha,q_delta}
parameters.V2=V2/epsilon;

angstore={};
spinstore={};

while 1
theta=rand(1,3)*2*pi;
phi=rand(1,3)*pi;
% theta=[0,0,0];
% phi=[1,0,0]*pi;

ansatz={[sin(phi(1))*cos(theta(1)),sin(phi(1))*sin(theta(1)),cos(phi(1))],...
    [sin(phi(2))*cos(theta(2)),sin(phi(2))*sin(theta(2)),cos(phi(2))],...
    [sin(phi(3))*cos(theta(3)),sin(phi(3))*sin(theta(3)),cos(phi(3))]};

Nsite=randi(3);
ansatz=ansatz(1:Nsite);
% param=mainTMD_spin('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',3,'nu',[1,2],'d',10,'spin_init',ansatz);
param.spin0={[0,0,1],ansatz{:}};
param.N=parameters.N;
param.energylist=parameters.energylist;
param.V1=parameters.V1;
param.V2=parameters.V2;

clear spinsav en gapsav
[energyall,wfall]=energyMF_init_2(param);
for i=1:999
[spin,gap,innergap]=spintexture(energyall,wfall,param);
[en(i),ave,V2deltaave]=totalenergy_2(energyall,wfall,param);
% fprintf("%d: gap:%0.8f meV E:%f meV innergap: %0.8f\n",i,1000*gap,1000*en(end),1000*innergap);
% disp([spin,angle(spin(:,2)+spin(:,3)*1i)*180/pi,angle(spin(:,4)+sqrt(spin(:,3).^2+spin(:,2).^2)*1i)*180/pi])
% % plot(en);
% ylabel('energy (eV)');
% drawnow;
% spinsav(:,:,i)=spin;
% % plot(squeeze(angle(spinsav(1,4,:)+sqrt(spinsav(1,3,:).^2+spinsav(1,2,:).^2)*1i)*180/pi));
% plot(en);
% gapsav(i)=gap;
if length(en)>1    
    if abs(en(end)-en(end-1))<1e-12
        break
    end
end
[energyall,wfall]=energyMF_2(ave,V2deltaave,param);
end
final=en(end);
ang=[0,0,0];
for ind=2:4
    ang(ind-1)=spin(ind,2:4)*spin(1,2:4)'/(norm(spin(ind,2:4))*norm(spin(1,2:4)));
end
exist=0;
for ind=1:length(angstore)
    if all(abs(ang-angstore{ind})<ones(1,3)*1e-2)
        exist=1;
        break;
    end
end

if exist==0
    disp([i,spin(:,1)',ang]);
    angstore={angstore{:},ang};
    spinstore={spinstore{:},spin};
end
end