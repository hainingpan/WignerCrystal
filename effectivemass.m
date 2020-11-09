function [nu,em]=effectivemass(theta)
parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'nu',[4,4],'d',60e-9*5.076e6,'Vz',0,'Ez',0);
%V:4.428
%w:20
n1=51;
n2=61;
h1=norm(parameters.bM1)/(n1);
h2=norm(parameters.bM2)/(n2);
% counter=1;
% clear kxlist kylist
kxlist=zeros(n1,n2);
kylist=zeros(n1,n2);
by=parameters.bM1-(parameters.bM1*parameters.bM2')*parameters.bM2/(norm(parameters.bM1)*norm(parameters.bM2));
bx=parameters.bM2;
for xindex=1:n1
    for yindex=1:n2
        u1=(2*xindex-n1-1)/(2*n1);
        u2=(2*yindex-n2-1)/(2*n2);
        klist=u1*parameters.bM1+u2*parameters.bM2;
        kxlist(xindex,yindex)=klist(1);
        kylist(xindex,yindex)=klist(2);
%         counter=counter+1;
    end
end

energyall=zeros(n1,n2);
parfor xindex=1:n1
    for yindex=1:n2
        val=energyTMD(kxlist(xindex,yindex),kylist(xindex,yindex),parameters);
        energyall(xindex,yindex)=val(1);
    end
end

e11=grad11(energyall,h1);
e22=grad22(energyall,h2);
e12=grad12(energyall,h1,h2);

U1x=[parameters.bM1'/norm(parameters.bM1),parameters.bM2'/norm(parameters.bM2)]; % k(x,y)=U1x k(1,2)
Ux1=inv(U1x); % k(1,2)=Ux1 k(x,y) see oneNote effective mass
jacobian=det(Ux1);
edet=(e11.*e22-e12.^2)*jacobian^2;

m0=0.511e6;

nu=nan*kxlist;
for i=1:n1
    for j=1:n2
        en=energyall(i,j);
        nu(i,j)=2*nnz(en<=energyall)/(n1*n2);        
    end
end
em=(edet(:)*m0^2);
fig=figure;
plot(nu(:),em,'.');
xlabel('\nu_{h}');
ylabel('m_x^{-1}*m_y^{-1}*m_e^2')
title(strcat('\theta=',num2str(theta),'{}^\circ'))
savefig(strcat('theta',num2str(theta),'.fig'));
saveas(fig,strcat('theta',num2str(theta),'.png'));
end