parameters=mainTMD('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',4,'nu',[1,3],'d',inf,'Vz',0);
tshell=3;
Ushell=35;
% [t,neighborlist]=t_calc_func(tshell,parameters);
U=U_calc_func_2(Ushell,parameters);


n=3;
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
hp=0;
tlist=-hp*[t{1:tshell+1}];
epsilon=1;
Ulist=real([U{1:Ushell+1}])/epsilon;

clear spinsav en
[energyall,wfall]=energyMF_init_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,parameters);
for i=1:200
[spin,gap]=spintexture(energyall,wfall,parameters);
en(i)=totalenergy_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,energyall,wfall,parameters);
fprintf("%d: gap:%0.8f meV E:%f meV\n",i,1000*gap,1000*en(end));
disp([spin,angle(spin(:,2)+spin(:,3)*1i)*180/pi,angle(spin(:,4)+sqrt(spin(:,3).^2+spin(:,2).^2)*1i)*180/pi])
plot(en);
ylabel('energy (eV)');
drawnow;
spinsav(:,:,i)=spin;
gapsav(i)=gap;
[energyall,wfall]=energyMF_2(kxlist,kylist,t_bond,tlist,U_bond,Ulist,energyall,wfall,parameters);
if length(en)>1    
    if abs(en(end)-en(end-1))<1e-8
        break
    end
end
end
final=en(end);
% save(sprintf('nu%d,%d_t%d_U%d_hp%d_ep%d.mat',parameters.nu(1),parameters.nu(2),tshell,Ushell,hp,epsilon),'en','spinsav','gapsav')
