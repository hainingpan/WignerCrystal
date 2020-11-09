function fermisurface(theta)
parameters=mainTMD_2('m',0.45,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',theta,'nu',[4,4],'d',60e-9*5.076e6,'Vz',0,'Ez',0);
%V:4.428
%w:20
n1=201;
n2=221;
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

kxlist=kxlist(:);
kylist=kylist(:);
energyall=energyall(:);
kxlist2=[];
kylist2=[];
energyall2=[];
for i=-1:1
    for j=-1:1
        kxlist2=[kxlist2;kxlist+i*parameters.bM1(1)+j*parameters.bM2(1)];
        kylist2=[kylist2;kylist+i*parameters.bM1(2)+j*parameters.bM2(2)];
        energyall2=[energyall2;energyall];
    end
end

rotate=@(x) [cos(x) -sin(x);sin(x) cos(x)]; %rotate anticlockwise
polygonpts=[(parameters.bM2-parameters.bM1)/3;
    (parameters.bM2-parameters.bM1)/3*rotate(pi/3);
    (parameters.bM2-parameters.bM1)/3*rotate(2*pi/3);
    (parameters.bM2-parameters.bM1)/3*rotate(3*pi/3);
    (parameters.bM2-parameters.bM1)/3*rotate(4*pi/3);
    (parameters.bM2-parameters.bM1)/3*rotate(5*pi/3);
    (parameters.bM2-parameters.bM1)/3];
[in,on]=inpolygon(kxlist2,kylist2,polygonpts(:,1),polygonpts(:,2));

kxlisthex=kxlist2(in);
kylisthex=kylist2(in);
energyallhex=energyall2(in);
[kxq,kyq]=meshgrid(linspace(min(kxlisthex),max(kxlisthex),100),linspace(min(kylisthex),max(kylisthex),100));
vq=griddata(kxlisthex,kylisthex,1000*energyallhex(:),kxq,kyq);
fig=figure;
mesh(kxq/norm(parameters.bM1),kyq/norm(parameters.bM1),vq);
hold on;
colorbar;
caxis(1000*[min(energyallhex) max(energyallhex)]);
[C,h]=contour(kxq/norm(parameters.bM1),kyq/norm(parameters.bM1),1000*vq-1000,15,'k');
h.ContourZLevel = -1000;
plot3(polygonpts(:,1)/norm(parameters.bM1),polygonpts(:,2)/norm(parameters.bM1),-1000.*ones(length(polygonpts(:,1))),'k');
view([0,-90]);
% view(3);
xlabel('k_x/|b_M|')
ylabel('k_y/|b_M|')
% zlabel('E(meV)')
daspect([1,1,40])
grid off;
title(strcat('E(meV) \theta=',num2str(theta),'{}^\circ'))
savefig(strcat('FS_theta',num2str(theta),'.fig'));
saveas(fig,strcat('FS_theta',num2str(theta),'.png'));




