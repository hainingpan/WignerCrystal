function parameters=mainTMD(varargin)
p=inputParser;
addParameter(p,'a',3.28e-10*5.076e6);
addParameter(p,'m',0.45);
addParameter(p,'theta',1.2);
addParameter(p,'V',8);
addParameter(p,'psi',-89.6);
addParameter(p,'w',-8.5);
addParameter(p,'Nmax',5);
addParameter(p,'Vz',0);
addParameter(p,'B',0);
addParameter(p,'nu',[1,1]);%filling factor number(#1) per site(#2); no spin degeneracy considered
parse(p,varargin{:});
parameters=struct('a',p.Results.a,'m',p.Results.m*0.511e6,'theta',p.Results.theta/360*2*pi,'V',p.Results.V*1e-3,'psi'...
    ,p.Results.psi/360*2*pi,'w',p.Results.w*1e-3,'Vz',p.Results.Vz*1e-3,'Nmax',p.Results.Nmax,'B',p.Results.B,'nu',p.Results.nu);
%Unit vectors
parameters.a1=parameters.a*[1,0];
parameters.a2=parameters.a*[1/2,sqrt(3)/2];
parameters.aM=parameters.a/parameters.theta;
parameters.aM1=parameters.aM*[0,-1];
parameters.aM2=parameters.aM*[sqrt(3)/2,-1/2];
% parameters.A1=parameters.aM*sqrt(3)*[1,0];
% parameters.A2=parameters.aM*sqrt(3)*[1/2,sqrt(3)/2];
%Reciprocal lattice
parameters.G1=[0,4*pi/(sqrt(3)*parameters.a)];
parameters.G2=[-2*pi/parameters.a,(2*pi)/(sqrt(3)*parameters.a)];
parameters.G3=[-2*pi/parameters.a,-(2*pi)/(sqrt(3)*parameters.a)];
parameters.G4=[0,-4*pi/(sqrt(3)*parameters.a)];
parameters.G5=[2*pi/parameters.a,-(2*pi)/(sqrt(3)*parameters.a)];
parameters.G6=[2*pi/parameters.a,(2*pi)/(sqrt(3)*parameters.a)];
parameters.b1=parameters.G5;    %reciprocal unit vector of small lattice 
parameters.b2=parameters.G1;    %reciprocal unit vector of small lattice
parameters.bM1=2*pi/parameters.aM*[-1/sqrt(3),-1];  %bM1=theta* cross(b1,z); reciprocal unit vector of Moire lattice
parameters.bM2=2*pi/parameters.aM*[2/sqrt(3),0];  %bM2=theta* cross(b2,z); reciprocal unit vector of Moire lattice
parameters.kp=2*pi/parameters.aM*[-1/sqrt(3),-1/3];  %k+=[-bM1.x/2,-bM1.x/2 tan(pi/6)]
parameters.kn=2*pi/parameters.aM*[-1/sqrt(3),1/3]; %k-=[-bM1.x/2,bM1.x/2 tan(pi/6)]
parameters.kpp=2*pi/parameters.aM*[1/sqrt(3),-1/3]; %k+'

%Potential term
Nrange=-parameters.Nmax:parameters.Nmax;
[h1index,h2index]=meshgrid(Nrange,Nrange);
parameters.h1index=h1index;
parameters.h2index=h2index;
[h1matX,h1matY]=meshgrid(h1index(:));
[h2matX,h2matY]=meshgrid(h2index(:));
h1mat=h1matX-h1matY;
h2mat=h2matX-h2matY;
parameters.DeltaTmat=reshape(arrayfun(@(h1,h2) DeltaT(h1,h2,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
parameters.DeltaTTmat=reshape(arrayfun(@(h1,h2) DeltaTT(h1,h2,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
parameters.Deltatmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,-1,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
parameters.Deltabmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,1,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);

%Zeeman field
mu_B=5.778e-5;
epsilon_z=(1+2+1/0.45)*mu_B;
Ez=epsilon_z*parameters.B;
parameters.Ez=Ez;

rotate=@(x) [cos(x) -sin(x);sin(x) cos(x)];
if parameters.nu==[1,2]
    ailist={[0,0],[-1,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,-1]};
    Qlist={[0,0],[1/2,0],[0,1/2],[1/2,1/2]};
    parameters.Q=cellfun(@(x) x(1)*parameters.bM1+x(2)*parameters.bM2,Qlist,'UniformOutput',0);
    parameters.Qindex=Qlist;
    parameters.Qindexmod=cellfun(@(x) [mod(x(1),1),mod(x(2),1)],Qlist,'UniformOutput',0);
    parameters.bm1=Qlist{2}(1)*parameters.bM1+Qlist{2}(2)*parameters.bM2;
    parameters.bm2=Qlist{3}(1)*parameters.bM1+Qlist{3}(2)*parameters.bM2;
end

if parameters.nu==[1,3]
    ailist={[0,0],[-2,1],[-4,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,0,0],[cos(2*pi/3),sin(2*pi/3),0],[cos(4*pi/3),sin(4*pi/3),0]};
    Qlist={[0,0],[1/3,0],[2/3,0],[0,1/3],[1/3,1/3],[2/3,1/3],[0,2/3],[1/3,2/3],[2/3,2/3]};
    parameters.Q=cellfun(@(x) x(1)*parameters.bM1+x(2)*parameters.bM2,Qlist,'UniformOutput',0);
    parameters.Qindex=Qlist;
    parameters.Qindexmod=cellfun(@(x) [mod(x(1),1),mod(x(2),1)],Qlist,'UniformOutput',0);
    parameters.bm1=Qlist{2}(1)*parameters.bM1+Qlist{2}(2)*parameters.bM2;
    parameters.bm2=Qlist{4}(1)*parameters.bM1+Qlist{4}(2)*parameters.bM2;
end

if parameters.nu==[2,3]
    ailist={[0,0],[-2,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,-1]};
    Qlist={[0,0],[-2/3,-1/3],[2/3,1/3]};
    parameters.Q=cellfun(@(x) x(1)*parameters.bM1+x(2)*parameters.bM2,Qlist,'UniformOutput',0);
    parameters.Qindex=Qlist;
    parameters.Qindexmod=cellfun(@(x) [mod(x(1),1),mod(x(2),1)],Qlist,'UniformOutput',0);
    parameters.bm1=Qlist{3}(1)*parameters.bM1+Qlist{3}(2)*parameters.bM2;
    parameters.bm2=parameters.bm1*rotate(pi/3);
end

if parameters.nu==[1,4]
    ailist={[0,0],[-2,2],[-4,4]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,0,0],[cos(2*pi/3),sin(2*pi/3),0],[cos(4*pi/3),sin(4*pi/3),0]};
    Qlist={[2/3,1/3],[1/3,1/6],[0,0],[-1/3,-1/6],[-2/3,-1/3],[-1,-1/2],...
        [2/3,1/3]+[5/6,4/6],[1/3,1/6]+[5/6,4/6],[5/6,4/6],[-1/3,-1/6]+[5/6,4/6],[-2/3,-1/3]+[5/6,4/6],[-1,-1/2]+[5/6,4/6]};
    parameters.Q=cellfun(@(x) x(1)*parameters.bM1+x(2)*parameters.bM2,Qlist,'UniformOutput',0);
    parameters.Qindex=Qlist;
    parameters.Qindexmod=cellfun(@(x) [mod(x(1),1),mod(x(2),1)],Qlist,'UniformOutput',0);
    parameters.bm1=Qlist{4}(1)*parameters.bM1+Qlist{4}(2)*parameters.bM2;
    parameters.bm2=parameters.bm1*rotate(pi/3);
end

if parameters.nu==[3,4]
    ailist={[0,0],[0,1],[-1,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[-1,0,0],[cos(-pi/3),sin(-pi/3),0],[cos(pi/3),sin(pi/3),0]};
    Qlist={[0,0],[1/2,0],[0,1/2],[1/2,1/2]};
    parameters.Q=cellfun(@(x) x(1)*parameters.bM1+x(2)*parameters.bM2,Qlist,'UniformOutput',0);
    parameters.Qindex=Qlist;
    parameters.Qindexmod=cellfun(@(x) [mod(x(1),1),mod(x(2),1)],Qlist,'UniformOutput',0);
    parameters.bm1=Qlist{2}(1)*parameters.bM1+Qlist{2}(2)*parameters.bM2;
    parameters.bm2=Qlist{3}(1)*parameters.bM1+Qlist{3}(2)*parameters.bM2;
end

if parameters.nu==[1,1]
    ailist={[0,0],[-1,1],[-2,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,0,0],[cos(2*pi/3),sin(2*pi/3),0],[cos(4*pi/3),sin(4*pi/3),0]};
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
    Qlist={[0,0],[-2/3,-1/3],[2/3,1/3]};
    parameters.Q=cellfun(@(x) x(1)*parameters.bM1+x(2)*parameters.bM2,Qlist,'UniformOutput',0);
    parameters.Qindex=Qlist;
    parameters.Qindexmod=cellfun(@(x) [mod(x(1),1),mod(x(2),1)],Qlist,'UniformOutput',0);
    parameters.bm1=Qlist{3}(1)*parameters.bM1+Qlist{3}(2)*parameters.bM2;
    parameters.bm2=parameters.bm1*rotate(pi/3);
end

end



