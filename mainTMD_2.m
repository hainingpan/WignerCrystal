function parameters=mainTMD_2(varargin)
p=inputParser;
addParameter(p,'a',3.28e-10*5.076e6);
addParameter(p,'m',0.45);
addParameter(p,'theta',1.2);
addParameter(p,'V',8);
addParameter(p,'psi',-89.6);
addParameter(p,'w',-8.5);
addParameter(p,'Nmax',5);
addParameter(p,'Vz',0);
addParameter(p,'Ez',0);
addParameter(p,'d',inf);
addParameter(p,'nu',[1,1]);%filling factor number(#1) per site(#2); no spin degeneracy considered
addParameter(p,'hole',1); %1: hole-like energy band ; -1: particle-like energy band
parse(p,varargin{:});
parameters=struct('a',p.Results.a,'m',p.Results.m*0.511e6,'theta',p.Results.theta/360*2*pi,'V',p.Results.V*1e-3,'psi'...
    ,p.Results.psi/360*2*pi,'w',p.Results.w*1e-3,'Vz',p.Results.Vz*1e-3,'Nmax',p.Results.Nmax,'Ez'...
    ,p.Results.Ez,'d',p.Results.d,'nu',p.Results.nu,'hole',p.Results.hole);
%Unit vectors
parameters.a1=parameters.a*[1,0];
parameters.a2=parameters.a*[1/2,sqrt(3)/2];
parameters.aM=parameters.a/parameters.theta;
% parameters.d=parameters.d*parameters.aM; %d in the unit of a moire
parameters.aM1=parameters.aM*[0,-1];
parameters.aM2=parameters.aM*[sqrt(3)/2,-1/2];
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
% mu_B=5.778e-5;
% epsilon_z=(1+2+1/0.45)*mu_B;
% Ez=epsilon_z*parameters.B;
% parameters.Ez=Ez;

rotate=@(x) [cos(x) -sin(x);sin(x) cos(x)]; %rotate anticlockwise
%For AF
if parameters.nu==[1,2] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,-1]};
%     parameters.spin0={[1,1,1]/sqrt(3),[1,-1,-1]/sqrt(3),[-1,1,-1]/sqrt(3),[-1,-1,1]/sqrt(3)};
   am1index=[0,2];
   am2index=[2,0];
end

%For F
if parameters.nu==[3,6] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,1]};
   am1index=[0,2];
   am2index=[2,0];
end

%For kagome-C2
if parameters.nu==[4,8] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1]*2/3,[0,0,1]*2/3,[0,0,1]*2/3};
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
    parameters.phi=0;
    parameters.s1=0.5;
   am1index=[0,2];
   am2index=[2,0];
end

%For kagome-T
if parameters.nu==[5,10] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1]*2/3,[0,0,1]*2/3,[0,0,1]*2/3};
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
    parameters.phi=pi/6;
    parameters.s1=1;
   am1index=[0,2];
   am2index=[2,0];
end

%For no Wigner Crystal
if parameters.nu==[6,12] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1]/2,[0,0,1]/2,[0,0,1]/2,[0,0,1]/2};

   am1index=[0,2];
   am2index=[2,0];
end

% % hexagonal 120 AF
% if parameters.nu==[6,12]
%     ailist={[0,1];[0,5];...
%             [-1,1];[-1,3];[-1,4];[-1,6];...
%             [-2,3];[-2,5];...
%             [-3,4];[-3,5];[-3,7];[-3,8];...
%             [-4,7];[-4,9];...
%             [-5,5];[-5,6];[-5,8];[-5,9];...
%             [0,0];[0,2];[0,3];[0,4];
%             [-1,2];[-1,5];...
%             [-2,2];[-2,4];[-2,6];[-2,7];...
%             [-3,3];[-3,6];...
%             [-4,4];[-4,5];[-4,6];[-4,8];...
%             [-5,7];[-5,10]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
% %     parameters.spin0={[1,-1,1]/sqrt(3),[1,1,-1]/sqrt(3),[-1,1,-1]/sqrt(3),[-1,1,1]/sqrt(3),[-1,-1,1]/sqrt(3),[1,-1,1]/sqrt(3)};
%     spin0=[1,0,0];spin1=[cos(pi/3),sin(pi/3),0];spin2=[cos(2*pi/3),sin(2*pi/3),0];
%     spin3=[cos(pi),sin(pi),0];spin4=[[cos(4*pi/3),sin(4*pi/3),0]];spin5=[cos(5*pi/3),sin(5*pi/3),0];
%     parameters.spin0={spin2,spin5,...
%                     spin5,spin0,spin3,spin2,...
%                     spin3,spin0,...
%                     spin0,spin3,spin4,spin1,...
%                     spin1,spin4,...
%                     spin2,spin5,spin4,spin1};
%     am1index=[0,6];
%     am2index=[-6,6];
% end

% % hexagonal F
% if parameters.nu==[12,24]
%     ailist={[0,0],[-2,2],[-3,3],[-5,5],[-2,1],[-4,3],...
%             [-1,1],[-4,4],[-3,2],[-1,2],[-2,3],[-3,4]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]};    
%     am1index=[-2,4];
%     am2index=[-4,2];
% end
% 
% if parameters.nu==[8,16]
%     ailist={[0,0],[-1,0],[-2,0],[-3,0],[-1,2],[-2,2],[-3,2],[-4,2],[-1,1],[-2,1],[-3,1],[-4,1],[-2,3],[-3,3],[-4,3],[-5,3]};    
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[-1,-1,1]/sqrt(3),[-1,-1,-1]/sqrt(3),[-1,1,-1]/sqrt(3),[-1,1,1]/sqrt(3),[1,-1,1]/sqrt(3),[1,-1,-1]/sqrt(3),[1,1,-1]/sqrt(3),[1,1,1]/sqrt(3)};
%    am1index=[-2,4];
%    am2index=[-4,0];
% end

if parameters.nu==[1,3]
    ailist={[0,0],[-2,1],[-4,2],[-1,1],[-2,2],[-3,2],[-3,1],[-2,0],[-1,0]};    
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,0,0],[cos(2*pi/3),sin(2*pi/3),0],[cos(4*pi/3),sin(4*pi/3),0]};
   am1index=[-3,3];
   am2index=[-3,0];
end

if parameters.nu==[2,6]
    ailist={[0,0],[-2,1],[-4,2],[-1,1],[-2,2],[-3,2],[-3,1],[-2,0],[-1,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[[0,0,-1]]};
   am1index=[-3,3];
   am2index=[-3,0];
end

%FM with 3 unit cells
if parameters.nu==[3,9]
    ailist={[0,0],[-2,2],[-1,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1]};
   am1index=[-1,2];
   am2index=[-2,1];
end

%Half honeycomb
if parameters.nu==[4,12]
    ailist={[0,0],[-2,2],[-1,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   am1index=[-1,2];
   am2index=[-2,1];
end


% %For tetrahedron spin
% if parameters.nu==[4,12]
%     ailist={[0,0],[-1,2],[-2,1],[-3,3],[-2,2],[-4,4],[-2,3],[-3,4],[-1,1],[-5,5],[-3,2],[-4,3]}; 
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[1,1,1]/sqrt(3),[1,-1,-1]/sqrt(3),[-1,1,-1]/sqrt(3),[-1,-1,1]/sqrt(3)};
%    am1index=[-2,4];
%    am2index=[-4,2];
% end

if parameters.nu==[2,3]
    ailist={[0,0],[-2,2],[-1,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,-1]};
   am1index=[-1,2];
   am2index=[-2,1];
end
% For 2/3 F
if parameters.nu==[6,9]
    ailist={[0,0],[-2,2],[-1,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,1]};
   am1index=[-1,2];
   am2index=[-2,1];
end
% For 2/3 120 AF
if parameters.nu==[8,12]
    ailist={[0,0],[-1,1],[-2,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,0,0],[cos(-2*pi/3),sin(-2*pi/3),0],[cos(-4*pi/3),sin(-4*pi/3),0]};
   am1index=[-1,2];
   am2index=[-2,1];
end

%
if parameters.nu==[4,6]
    ailist={[0,0],[-2,2],[-1,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);

   am1index=[-1,2];
   am2index=[-2,1];
end
% AF 120
if parameters.nu==[1,4]
    ailist={[0,0],[-2,2],[-4,4],[-1,2],[-2,3],[-3,4],[-1,1],[-3,3],[-5,5],[-2,1],[-3,2],[-4,3]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,0,0],[cos(2*pi/3),sin(2*pi/3),0],[cos(4*pi/3),sin(4*pi/3),0]};
   am1index=[-2,4];
   am2index=[-4,2];
end
% For ferromagnetic
if parameters.nu==[2,8]
    ailist={[0,0],[-2,2],[-4,4],[-1,2],[-2,3],[-3,4],[-1,1],[-3,3],[-5,5],[-2,1],[-3,2],[-4,3]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
   am1index=[-2,4];
   am2index=[-4,2];
end

%  %canted AF as ansatz
% if parameters.nu==[3,12]
%     ailist={[0,0],[-2,2],[-4,4],[-1,2],[-2,3],[-3,4],[-1,1],[-3,3],[-5,5],[-2,1],[-3,2],[-4,3]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[1,0,cos(pi/4)],[cos(2*pi/3),sin(2*pi/3),cos(pi/4)],[cos(4*pi/3),sin(4*pi/3),cos(pi/4)]};
%    am1index=[-2,4];
%    am2index=[-4,2];
% end

% For kagome lattice AF
if parameters.nu==[3,12]
    ailist={[0,0],[-1,2],[-2,1],[-1,1],[-2,2],[-3,3],[-4,4],[-5,5],[-3,2],[-4,3],[-2,3],[-3,4]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   parameters.spin0={[cos(-5*pi/6),sin(-5*pi/6),0],[cos(-pi/6),sin(-pi/6),0],[cos(pi/2),sin(pi/2),0]};
   am1index=[-2,4];
   am2index=[-4,2];
end

% For kagome lattice F
if parameters.nu==[4,16]
    ailist={[0,0],[-1,2],[-2,1],[-1,1],[-2,2],[-3,3],[-4,4],[-5,5],[-3,2],[-4,3],[-2,3],[-3,4]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
   am1index=[-2,4];
   am2index=[-4,2];
end

% For stripe AF
if parameters.nu==[5,20]
    ailist={[0,0];[-1,2];[-1,1];[-2,2];[-3,3];[-2,3];[-3,4];[-4,5]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   parameters.spin0={[0,0,1],[0,0,-1]};
   am1index=[-2,4];
   am2index=[-3,2];
end

% % For stripe F
if parameters.nu==[6,24]
    ailist={[0,0];[-1,2];[-1,1];[-2,2];[-3,3];[-2,3];[-3,4];[-4,5]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   parameters.spin0={[0,0,1],[0,0,1]};
   am1index=[-2,4];
   am2index=[-3,2];
end

%For kagome-C2
if parameters.nu==[7,28] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.phi=0;
    parameters.s1=0.5;
   am1index=[0,2];
   am2index=[2,0];
end

%For kagome-T
if parameters.nu==[8,32] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.phi=pi/6;
    parameters.s1=1;
   am1index=[0,2];
   am2index=[2,0];
end

%For tetrahedron spin
if parameters.nu==[2,4] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,1,1]/sqrt(3),[1,-1,-1]/sqrt(3),[-1,1,-1]/sqrt(3),[-1,-1,1]/sqrt(3)};
%     parameters.spin0={[1,1,1]/sqrt(3),[1,-1,-1]/sqrt(3),-[-1,1,-1]/sqrt(3),-[-1,-1,1]/sqrt(3)}; % two in two out
%     parameters.spin0={[1,1,1]/sqrt(3),[1,-1,-1]/sqrt(3),[-1,1,-1]/sqrt(3),-[-1,-1,1]/sqrt(3)}; % three in one out

%     parameters.spin0={[0,0,1];[0,0,-1]/3;[0,0,-1]/3;[0,0,-1]/3};
   am1index=[0,2];
   am2index=[2,0];
end

%Kagome AF
if parameters.nu==[3,4]
    ailist={[0,0],[0,1],[-1,1],[-1,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[-1,0,0],[cos(-pi/3),sin(-pi/3),0],[cos(pi/3),sin(pi/3),0]};
   am1index=[0,2];
   am2index=[-2,2];
end

%Kagome F
if parameters.nu==[6,8]
    ailist={[0,0],[0,1],[-1,1],[-1,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
   am1index=[0,2];
   am2index=[-2,2];
end

% For complimentary-kagome lattice AF
if parameters.nu==[9,12]
    ailist={[-1,1],[-2,2],[-3,3],[-4,4],[-5,5],[-3,2],[-4,3],[-2,3],[-3,4],...
        [0,0],[-1,2],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    spinup=[0,0,1];spindown=[0,0,-1];
    parameters.spin0={spinup,spindown,spinup,spindown,spinup,spindown,spindown,spindown,spindown,spindown,spindown,spindown};
    am1index=[-2,4];
   am2index=[-4,2];
end


% For complimentary-kagome lattice F
if parameters.nu==[12,16]
    ailist={[-1,1],[-2,2],[-3,3],[-4,4],[-5,5],[-3,2],[-4,3],[-2,3],[-3,4],...
        [0,0],[-1,2],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   parameters.spin0={[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]};
   am1index=[-2,4];
   am2index=[-4,2];
end

% For complimentary-stripe AF
% if parameters.nu==[15,20]
%     ailist={[-1,0];[-3,0];[-5,0];[-1,1];[-2,1];[-3,1];[-4,1];[-5,1];[-6,1];...
%             [0,0];[-2,0];[-4,0]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     spin0=[1,0,0];spin1=[cos(2*pi/3),sin(2*pi/3),0];spin2=[cos(4*pi/3),sin(4*pi/3),0];
%     spinfunc=@(x) [cos(x*pi/180),sin(x*pi/180),0];
%     spinup=[0,0,1];spindown=[0,0,-1];
% %    parameters.spin0={spin2,spin0,spin1,spin1,spin0,spin2,spin1,spin0,spin2};
% %     x=rand(1,3)*360;
% %     x=[120,0,240];
% %     parameters.spin0={spinfunc(0),spinfunc(180),spinfunc(0),spinfunc(120),spinfunc(240),spinfunc(60),spinfunc(-60),spinfunc(-240),spinfunc(-120)};
% %     parameters.spin0={spinup,spinup,spinup,spindown,spindown,spindown,spindown,spindown,spindown};
%    am1index=[-1,2];
%    am2index=[-6,0];
% end
 
%Skew
% if parameters.nu==[15,20]
%     ailist={[-1,1],[-2,2],[-3,3],[-2,3],[-3,4],[-4,5],[0,0],[-1,2]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     spinfunc=@(x) [cos(x*pi/180),sin(x*pi/180),0];
%     parameters.spin0={spinfunc(-90),spinfunc(180),spinfunc(-90),spinfunc(90),spinfunc(0),spinfunc(90)};
%    am1index=[-2,4];
%    am2index=[-3,2];
% end
%rect
if parameters.nu==[15,20]
    ailist={[0,0],[-2,0],[-1,1],[-2,1],[-3,1],[-4,1],[-1,0],[-3,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    spinfunc=@(x) [cos(x*pi/180),sin(x*pi/180),0];
   parameters.spin0={spinfunc(0),spinfunc(180),spinfunc(240),spinfunc(60),spinfunc(-60),spinfunc(-240)};
   am1index=[-1,2];
   am2index=[-4,0];
end


% if parameters.nu==[15,20]
%     ailist={[0,0],[-2,0],[-1,1],[-2,1],[-3,1],[-4,1],[-1,2],[-3,2],[-2,3],[-3,3],[-4,3],[-5,3],...
%         [-1,0],[-3,0],[-2,2],[-4,2]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     spin1=[1,1,1]/sqrt(3);spin2=[1,-1,-1]/sqrt(3);spin3=[-1,1,-1]/sqrt(3);spin4=[-1,-1,1]/sqrt(3);
%    parameters.spin0={spin1,-spin1,spin2,-spin2,-spin4,spin4,spin3,-spin3,spin4,-spin4,-spin2,spin2};
%    am1index=[-4,0];
%    am2index=[-2,4];
% end

% For complimentary-stripe F
% if parameters.nu==[18,24]
%     ailist={[-1,0];[-3,0];[-5,0];[-1,1];[-2,1];[-3,1];[-4,1];[-5,1];[-6,1];...
%             [0,0];[-2,0];[-4,0]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%    parameters.spin0={[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]};
%    am1index=[-1,2];
%    am2index=[-6,0];
% end
if parameters.nu==[18,24]
    ailist={[-1,1],[-2,1],[-2,2],[0,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
   am1index=[-1,2];
   am2index=[-2,0];
end

%For kagome-C2
if parameters.nu==[21,28]
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.phi=0;
    parameters.s1=0.5;
    am1index=[0,2];
   am2index=[-2,2];
end

%For kagome-T
if parameters.nu==[24,32] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.phi=pi/6;
    parameters.s1=1;
    am1index=[0,2];
   am2index=[-2,2];
end

%For tetrahedron spin
if parameters.nu==[4,4] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,1,1]/sqrt(3),[1,-1,-1]/sqrt(3),[-1,1,-1]/sqrt(3),[-1,-1,1]/sqrt(3)};
   am1index=[0,2];
   am2index=[2,0];
end

%For Wigner crystal
if parameters.nu==[1,5] 
    ailist={[0,0];[-1,1];[-2,1];[-2,2];[-3,2]};
%     ailist={[0,0];[-1,0];[-2,0];[-3,0];[-4,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1]};
%     parameters.spin0={[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]};
    am1index=[-1,2];
    am2index=[-3,1];
%     am1index=[-5,0];
%     am2index=[-1,1];
end

%120 AF
if parameters.nu==[2,10] 
    ailist={[0,0];[-1,2];[-2,4];...
            [-1,1];[0,1];[-2,2];[0,2];[-3,3];[-2,3];[-1,3];[0,3];[-3,4];[-1,4];[-3,5];[-2,5]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,0,0];[cos(2*pi/3),sin(2*pi/3),0];[cos(4*pi/3),sin(4*pi/3),0]};
    am1index=[-4,3];
    am2index=[1,3];
end

%180 AF-skewed
if parameters.nu==[3,15] 
    ailist={[0,0];[-3,1];[-1,2];[-4,3];...
        [-1,1];[-2,1];[-2,2];[-3,2];[-4,2];[-5,2];[-2,3];[-3,3];[-5,3];[-6,3];[-3,4];[-4,4];[-5,4];[-6,4];[-6,5];[-7,5]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,-1],[0,0,-1],[0,0,1]};
    am1index=[-2,4];
    am2index=[-6,2];
end


%180 AF-rectangular
if parameters.nu==[4,20] 
    ailist={[0,0];[-3,1];[-1,2];[-4,3];...
        [-1,0];[-2,0];[-3,0];[-4,0];[-1,1];[-2,1];[-4,1];[-5,1];[-2,2];[-3,2];[-4,2];[-5,2];[-2,3];[-3,3];[-5,3];[-6,3]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,-1],[0,0,-1],[0,0,1]};
    am1index=[-2,4];
    am2index=[-5,0];
end



% For spin texture
if parameters.nu==[8,8] 
    ailist={[0,0];[0,1];[0,2];[0,3];[1,1];[1,2];[-1,2];[-1,3]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,-1];[1,0,0];[0,0,1];[1,0,0];[cos(4*pi/3),sin(4*pi/3),0];...
%         [cos(2*pi/3),sin(2*pi/3),0];[cos(2*pi/3),sin(2*pi/3),0];[cos(4*pi/3),sin(4*pi/3),0]};
    spin1=[0,0,1];spin2=[0,0,-1];
    spin3=[1,0,0];spin4=[cos(2*pi/3),sin(2*pi/3),0];spin5=[cos(4*pi/3),sin(4*pi/3),0];
    parameters.spin0={[0,0,-1];[1,0,0];[0,0,1];[-1,0,0];[cos(5*pi/3),sin(5*pi/3),0];...
        [cos(4*pi/3),sin(4*pi/3),0];[cos(pi/3),sin(pi/3),0];[cos(2*pi/3),sin(2*pi/3),0]};
%     parameters.spin0={spin3;spin2;spin3;spin1;spin5;spin5;spin4;spin4};
    am1index=[2,1];
    am2index=[-2,3];
end
% 
% %Wigner crystal-AF
% if parameters.nu==[2,5] 
%     ailist={[-1,1];[-2,1];[-2,2];[-3,2];[0,0]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1],[0,0,-1]};
%     am1index=[-1,2];
%     am2index=[-3,1];
% end

% %Wigner crystal-F
% if parameters.nu==[2,5] 
%     ailist={[-1,1];[-2,1];[-2,2];[-3,2];[0,0]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1],[0,0,1]};
%     am1index=[-1,2];
%     am2index=[-3,1];
% end

% 180 AF for 10 unit cell
if parameters.nu==[2,5] 
    ailist={[0,0];[-1,2];[-2,1];[-3,3];...
            [-1,1];[-2,3];[-2,2];[-3,4];[-3,2];[-4,4]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,-1],[0,0,-1],[0,0,1]};
    am1index=[-2,4];
    am2index=[-3,1];
end

% 120 AF for 30 unit cell
% if parameters.nu==[2,5] 
%     ailist={[0,0];[4,0];[1,1];[2,1];[-1,2];[3,2];[0,3];[1,3];[-2,4];[2,4];[-1,5];[0,5];...
%             [1,0];[2,0];[3,0];[0,1];[3,1];[4,1];[0,2];[1,2];[2,2];[-1,3];[2,3];[3,3];[-1,4];[0,4];[1,4];[-2,5];[1,5];[2,5]};
%     spin0=[1,0,0];spin1=[cos(pi/3),sin(pi/3),0];spin2=[cos(2*pi/3),sin(2*pi/3),0];
%     spin3=[cos(pi),sin(pi),0];spin4=[[cos(4*pi/3),sin(4*pi/3),0]];spin5=[cos(5*pi/3),sin(5*pi/3),0];
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={spin0,spin3,spin4,spin1,spin2,spin5,spin0,spin3,spin4,spin1,spin2,spin5};
%     am1index=[5,0];
%     am2index=[-3,6];
% end

 
% if parameters.nu==[3,5] 
% %     ailist={[0,0];[-1,1];[-2,1];[-2,2];[-3,2]};
%     ailist={[-3,2];[-1,1];[-2,1];[-2,2];[0,0]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
%     am1index=[-1,2];
%     am2index=[-3,1];
% end


%For spin texture
if parameters.nu==[6,10] 
    ailist={[0,0];[-2,2];[-3,2];[-3,1];[-5,3];[-6,3];[-1,1];[-2,1];[-4,2];[-5,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,-1],[0,0,1],[0,0,-1],[0,0,1],[0,0,-1]};
    am1index=[-1,2];
    am2index=[-6,2];
end

%For Wigner Crystal
if parameters.nu==[4,5] 
    ailist={[-1,1];[-2,1];[-2,2];[-3,2];[0,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,1],[0,0,1],[0,0,1]};
    am1index=[-1,2];
    am2index=[-3,1];
end

%120 AF for 10 unit cell
if parameters.nu==[8,10] 
    ailist={[0,0];[2,0];[3,0];[4,0];[0,1];[1,1];[2,1];[4,1];[1,0];[3,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    spin0=[1,0,0];spin1=[cos(pi/3),sin(pi/3),0];spin2=[cos(2*pi/3),sin(2*pi/3),0];
    spin3=[cos(pi),sin(pi),0];spin4=[[cos(4*pi/3),sin(4*pi/3),0]];spin5=[cos(5*pi/3),sin(5*pi/3),0];
    parameters.spin0={spin2,spin5,spin3,spin0,spin0,spin3,spin1,spin4};
%     spin=@(x) [cos(x),sin(x),0];
%     theta1=119.8143/180*pi;
%     theta2=179.7115/180*pi;
%     parameters.spin0={spin(theta1),-spin(theta1),spin(theta2),-spin(theta2),spin(-theta2),-spin(-theta2),-spin(-theta1),spin(-theta1)};
    am1index=[-1,2];
    am2index=[5,0];
end

% % 120 AF for 15 unit cell
if parameters.nu==[12,15] 
    ailist={[0,0];[0,1];[1,1];[2,1];[4,1];[5,1];[6,1];...
            [1,2];[2,2];[3,2];[4,2];[6,2];...
            [3,1];[0,2];[5,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    spin0=[1,0,0];spin1=[cos(2*pi/3),sin(2*pi/3),0];spin2=[cos(4*pi/3),sin(4*pi/3),0];
    parameters.spin0={spin0,spin2,spin0,spin1,spin0,spin1,spin2,...
                    spin2,spin0,spin1,spin2,spin1};
    am1index=[7,1];
    am2index=[-1,2];
end

%No Wigner Crystal
if parameters.nu==[16,20] 
    ailist={[-1,1];[-2,1];[-2,2];[-3,2];[0,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1]*4/5,[0,0,1]*4/5,[0,0,1]*4/5,[0,0,1]*4/5,[0,0,1]*4/5};
    am1index=[-1,2];
    am2index=[-3,1];
end


if parameters.nu==[1,6] 
    ailist={[0,0];[1,0];[2,0];[0,1];[1,1];[2,1]};
%     ailist={[0,0];[1,0];[2,0];[1,1];[2,1];[3,1]};
%     ailist={[0,0];[0,1];[1,1];[2,1];[1,2];[2,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1]};
    am1index=[3,0];
    am2index=[0,2];
    
%     am1index=[3,0];
%     am2index=[1,2];
    
%     am1index=[3,1];
%     am2index=[0,2];
end

% if parameters.nu==[2,6] 
% %     ailist={[0,0];[1,0];[2,0];[0,1];[1,1];[2,1]};
%     ailist={[1,1];[2,1];[0,0];[1,0];[2,0];[3,1]};
% %     ailist={[0,0];[1,1];[0,1];[2,1];[1,2];[2,2]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1],[[0,0,-1]]};
% %     am1index=[3,0];
% %     am2index=[0,2];
%     
%     am1index=[3,0];
%     am2index=[1,2];
%     
% %     am1index=[3,1];
% %     am2index=[0,2];
% end

% if parameters.nu==[3,6] 
% %     ailist={[0,0];[1,0];[2,0];[0,1];[1,1];[2,1]};
% %     ailist={[0,0];[1,0];[2,0];[1,1];[2,1];[3,1]};
%     ailist={[0,1];[2,1];[1,2];[0,0];[1,1];[2,2]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
% %     am1index=[3,0];
% %     am2index=[0,2];
%     
% %     am1index=[3,0];
% %     am2index=[1,2];
%     
%     am1index=[3,1];
%     am2index=[0,2];
% end

% if parameters.nu==[4,6] 
% %     ailist={[1,0];[0,1];[1,1];[2,1];[0,0];[2,0]};
% %     ailist={[1,0];[2,0];[1,1];[2,1];[3,1];[0,0]};
%     ailist={[0,0];[0,1];[1,1];[2,1];[1,2];[2,2]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1],[0,0,1]};
% %     am1index=[3,0];
% %     am2index=[0,2];
%     
% %     am1index=[3,0];
% %     am2index=[1,2];
%     
%     am1index=[3,1];
%     am2index=[0,2];
% end

% FM Insulator
if parameters.nu==[5,6] 
    ailist={[0,0];[1,0];[0,1];[1,1];[2,1];[2,0]};
%     ailist={[1,0];[2,0];[1,1];[2,1];[3,1];[0,0]};
%     ailist={[0,0];[0,1];[1,1];[2,1];[1,2];[2,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]};
    am1index=[3,0];
    am2index=[0,2];
    
%     am1index=[3,0];
%     am2index=[1,2];
    
%     am1index=[3,1];
%     am2index=[0,2];
end

%FM insulator
if parameters.nu==[7,6] 
    ailist={[0,0];[1,0];[0,1];[1,1];[2,1];[2,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]};
    am1index=[3,0];
    am2index=[0,2];
end

% FM metalic
if parameters.nu==[5,6]*2 
    ailist={[0,0];[1,0];[0,1];[1,1];[2,1];[2,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1]*5/6,[0,0,1]*5/6,[0,0,1]*5/6,[0,0,1]*5/6,[0,0,1]*5/6,[0,0,1]*5/6};
    parameters.spin0={[0,0,1]*0,[0,0,1]*0,[0,0,1]*(0),[0,0,1]*(0),[0,0,1]*0,[0,0,1]*(0)};
    am1index=[3,0];
    am2index=[0,2];
end

% FM metalic
if parameters.nu==[7,6]*2 
    ailist={[0,0];[1,0];[0,1];[1,1];[2,1];[2,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1]*5/6,[0,0,1]*5/6,[0,0,1]*5/6,[0,0,1]*5/6,[0,0,1]*5/6,[0,0,1]*5/6};
    am1index=[3,0];
    am2index=[0,2];
end
% AF (use the typical in-plane 120 AF)
if parameters.nu==[5,6]*3 
    ailist={[0,0];[0,1];[0,2];[0,3];[0,4];[0,5];
            [-1,1];[-1,3];[-1,5];
            [-2,2];[-2,3];[-2,4];[-2,5];[-2,6];[-2,7];
            [-1,2];[-1,4];[-1,6]};
    
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    s0=[1,0,0];
    s1=[cos(2*pi/3),sin(2*pi/3),0];
    s2=[cos(-2*pi/3),sin(-2*pi/3),0];
    parameters.spin0={s1;s0;s2;s1;s0;s2;
                    s2;s0;s1;
                    s0;s2;s1;s0;s2;s1};
    am1index=[-3,3];
    am2index=[0,6];
end

%normal metal
if parameters.nu==[5,6]*4
    ailist={[0,0];[1,0];[0,1];[1,1];[2,1];[2,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]};
    am1index=[3,0];
    am2index=[0,2];

end


% % For Wigner Crystal
% if parameters.nu==[1,7] 
%     ailist={[0,0];[-1,1];[-1,2];[-2,2];[-2,3];[-3,3];[-3,4]};
% %     ailist={[0,0];[1,1];[2,1];[3,1];[2,2];[3,2];[4,2]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1]};
%     am1index=[-1,3];
%     am2index=[-3,2];
% %     am1index=[4,1];
% %     am2index=[1,2];
% end

% For spin texture
if parameters.nu==[1,7] 
    ailist={[0,0];[-1,3];[-2,6];[0,1];[0,2];[0,3];[0,4];[0,5];[-1,2];[-1,4];[-1,5];[-1,6];[-2,3];[-2,4];[-2,5];[-2,7];...
        [-3,4];[-3,5];[-3,6];[-3,7];[-3,8]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[1,0,0],[cos(2*pi/3),sin(2*pi/3),0],[cos(4*pi/3),sin(4*pi/3),0]};
    am1index=[1,4];
    am2index=[-4,5];
end

if parameters.nu==[2,7] 
    ailist={[-2,2];[-2,3];[0,0];[-1,1];[-1,2];[-3,3];[-3,4]};
%     ailist={[1,1];[2,1];[3,1];[2,2];[3,2];[4,2];[0,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,1]};
    am1index=[-1,3];
    am2index=[-3,2];
%     am1index=[4,1];
%     am2index=[1,2];
end

% For Wigner Crystal
if parameters.nu==[6,7] 
    ailist={[-1,1];[-1,2];[-2,2];[-2,3];[-3,3];[-3,4];[0,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1]};
    am1index=[-1,3];
    am2index=[-3,2];
%     am1index=[4,1];
%     am2index=[1,2];
end


% For spin texture
% if parameters.nu==[6,7] 
%     ailist={[0,1];[0,2];[0,3];[0,4];[0,5];[-1,2];[-1,4];[-1,5];[-1,6];[-2,3];[-2,4];[-2,5];[-2,7];...
%         [-3,4];[-3,5];[-3,6];[-3,7];[-3,8];[0,0];[-1,3];[-2,6]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     spin0=[1,0,0];
%     spin1=[cos(2*pi/3),sin(2*pi/3),0];
%     spin2=[cos(4*pi/3),sin(4*pi/3),0];
%     parameters.spin0={spin2,spin1,spin0,spin2,spin1,...
%                        spin0,spin1,spin0,spin2,...,
%                         spin1,spin0,spin2,spin0,...
%                         spin2,spin1,spin0,spin2,spin1};
%     am1index=[1,4];
%     am2index=[-4,5];
% end
% 
% % For spin texture
% if parameters.nu==[6,7] 
%     ailist={[0,0];[0,1];[0,2];[0,4];[0,5];[0,6];...
%             [1,1];[1,2];[1,3];[1,4];[1,5];[1,6];...
%             [2,1];[2,2];[2,3];[2,5];[2,6];[2,7];...
%             [3,2];[3,3];[3,4];...
%             [-1,3];[-1,4];[-1,5];...
%             [0,3];[1,7];[2,4];[-1,6]};
%     parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
% 
%     spin1=[1,1,1]/sqrt(3);
%     spin2=[1,-1,-1]/sqrt(3);
%     spin3=[-1,1,-1]/sqrt(3);
%     spin4=[-1,-1,1]/sqrt(3);
%     parameters.spin0={spin1,spin4,spin1,spin1,spin4,spin1,...
%                        spin2,spin3,spin2,spin3,spin2,spin3,...
%                         spin4,spin1,spin4,spin4,spin1,spin4,...
%                         spin3,spin2,spin3,...
%                         spin2,spin3,spin2};
%     am1index=[4,2];
%     am2index=[-2,6];
% end

%FM metallic
if parameters.nu==[12,14] 
    ailist={[-1,1];[-1,2];[-2,2];[-2,3];[-3,3];[-3,4];[0,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1]*6/7;[0,0,1]*6/7;[0,0,1]*6/7;[0,0,1]*6/7;[0,0,1]*6/7;[0,0,1]*6/7;[0,0,1]*6/7};
    am1index=[-1,3];
    am2index=[-3,2];
%     am1index=[4,1];
%     am2index=[1,2];
end

%AF
if parameters.nu==[5,12] 
    ailist={[0,0];[-1,0];[1,-1];[0,1];[-2,2];
        [0,-1];[1,0];[-1,1];[1,1];[0,2];[-1,2];[-2,1]};

    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);

    parameters.spin0={[0,0,0];[cos(pi/2),sin(pi/2),0];[cos(7/6*pi),sin(7/6*pi),0];[cos(-1/6*pi),sin(-1/6*pi),0];[cos(7/6*pi),sin(7/6*pi),0]};
    am1index=[2,2];
    am2index=[-2,4];
end

% FM
if parameters.nu==[10,24] 
    ailist={[0,0];[-1,0];[1,-1];[0,1];[-2,2];...
        [0,-1];[1,0];[-1,1];[1,1];[0,2];[-1,2];[-2,1]};

    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);

    parameters.spin0={[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1]};
    am1index=[2,2];
    am2index=[-2,4];
end

% no WC
if parameters.nu==[15,36] 
    ailist={[0,0];[-1,0];[1,-1];[0,1];[-2,2];...
        [0,-1];[1,0];[-1,1];[1,1];[0,2];[-1,2];[-2,1]};

    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);

    parameters.spin0={[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12;[0,0,1]*5/12};
    am1index=[2,2];
    am2index=[-2,4];
end

%AF -- with sqrt(3)* sqrt(3) triangular
if parameters.nu==[20,48] 
    occupied={[0,0];[-1,0];[1,-1];[0,1];[-2,2]};
    occupied1=cellfun(@(x) x+[2,2],occupied,'UniformOutput',false);
    occupied2=cellfun(@(x) x+[-2,4],occupied,'UniformOutput',false);
    unoccupied={[0,-1];[1,0];[-1,1];[1,1];[0,2];[-1,2];[-2,1]};
    unccoupied1=cellfun(@(x) x+[2,2],unoccupied,'UniformOutput',false);
    unccoupied2=cellfun(@(x) x+[-2,4],unoccupied,'UniformOutput',false);

    ailist={occupied{:},occupied1{:},occupied2{:},unoccupied{:},unccoupied1{:},unccoupied2{:}};

    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);

    parameters.spin0=...
        {[0,0,0];[cos(pi/2),sin(pi/2),0];[cos(-1/6*pi),sin(-1/6*pi),0];[cos(-5/6*pi),sin(-5/6*pi),0];[cos(-1/6*pi),sin(-1/6*pi),0];
        [0,0,0];[cos(-1/6*pi),sin(-1/6*pi),0];[cos(-5/6*pi),sin(-5/6*pi),0];[cos(pi/2),sin(pi/2),0];[cos(-5/6*pi),sin(-5/6*pi),0];
        [0,0,0];[cos(-5/6*pi),sin(-5/6*pi),0];[cos(pi/2),sin(pi/2),0];[cos(-1/6*pi),sin(-1/6*pi),0];[cos(pi/2),sin(pi/2),0]};
    am1index=[0,6];
    am2index=[-6,6];
end
% FM rectangular- akin to 1/2
if parameters.nu==[5,12]*5
    ailist={[0,0];[1,-1];[-2,2];[-1,1];[1,1];[0,2];
        [-1,0];[0,1];[0,-1];[1,0];[-1,2];[-2,1]};

    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);

    parameters.spin0={[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1]};
    am1index=[2,2];
    am2index=[-2,4];
end

%FM kagome- akin to 1/2
if parameters.nu==[5,12]*6
%     ailist={[0,0];[-1,0];[1,-1];[0,1];[-2,2];
%         [0,-1];[1,0];[-1,1];[1,1];[-1,2];[-2,1]};
%      ailist={[-1,0];[1,-1];[0,1];
%         [0,-1];[1,0];[-1,1];[1,1];[0,2];[-1,2];[-2,1]};
    ailist={[-1,0];[1,-1];[0,1];[0,-1];[1,0];[-1,1];[1,1];[-1,2];[-2,1];
        [0,0];[-2,2];[0,2]};

    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);

    parameters.spin0={[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1];[0,0,1]};
    am1index=[2,2];
    am2index=[-2,4];
end

% FM
if parameters.nu==[6,5] 
    ailist={[-1,1];[-2,1];[-2,2];[-3,2];[0,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1],[0,0,1],[0,0,1],[0,0,1]};
    am1index=[-1,2];
    am2index=[-3,1];
end

% AF
if parameters.nu==[6,5]*2
    ailist={[-1,1];[-2,1];[-2,2];[-3,2];[0,0]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
    parameters.spin0={[0,0,1]*6/5,[0,0,1]*6/5,[0,0,1]*6/5,[0,0,1]*6/5,[0,0,1]*6/5};
    am1index=[-1,2];
    am2index=[-3,1];
end



%FM
if parameters.nu==[1,1]
    ailist={[0,0],[-1,1],[-2,2]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[1,0,0],[cos(-2*pi/3),sin(-2*pi/3),0],[cos(-4*pi/3),sin(-4*pi/3),0]};
    parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
   am1index=[-1,2];
   am2index=[-2,1];
end

%Kagome-C2
if parameters.nu==[2,2] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1]*2/3,[0,0,1]*2/3,[0,0,1]*2/3};
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
    parameters.phi=0;
    parameters.s1=0.5;
   am1index=[0,2];
   am2index=[2,0];
end

%Kagome-T
if parameters.nu==[4,4] 
    ailist={[0,0],[-1,0],[-1,1],[-2,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
%     parameters.spin0={[0,0,1]*2/3,[0,0,1]*2/3,[0,0,1]*2/3};
%     parameters.spin0={[0,0,1],[0,0,1],[0,0,1]};
    parameters.phi=pi/6;
    parameters.s1=1;
   am1index=[0,2];
   am2index=[2,0];
end

%Haldane
if parameters.nu==[6,6] 
    ailist={[0,0],[-2,2],[-1,1]};
    parameters.inner=cellfun(@(x) x(1)*parameters.aM1+x(2)*parameters.aM2,ailist,'UniformOutput',0);
   am1index=[-1,2];
   am2index=[-2,1];
end

am1=am1index*[parameters.aM1;parameters.aM2];
am2=am2index*[parameters.aM1;parameters.aM2];
parameters.bm1=2*pi/abs(length(ailist)*parameters.aM^2*sqrt(3)/2)*am1*rotate(-pi/2);
parameters.bm2=2*pi/abs(length(ailist)*parameters.aM^2*sqrt(3)/2)*am2*rotate(pi/2);
qindex=[parameters.bm1;parameters.bm2]/[parameters.bM1;parameters.bM2];
manhattandist=1;
Qlist={[0,0]};
done=0;
while 1  %Very ugly code
    for i=-manhattandist:manhattandist
        if abs(i)==manhattandist
            jlist=[0];
        else
            jlist=[-(manhattandist-abs(i)),(manhattandist-abs(i))];
        end
        for j=jlist
        qtmp=[i,j]*qindex;
        if all(qtmp>=-1e-10) & all(qtmp<1-1e-10)
            Qlist=[Qlist,qtmp];
        end
        if length(Qlist)==length(ailist)
            done=1;
            break;
        end
        end
    end
    if done==1
        break;
    end
    manhattandist=manhattandist+1;
end
parameters.Qindex=Qlist;
parameters.Q=cellfun(@(x) x(1)*parameters.bM1+x(2)*parameters.bM2,Qlist,'UniformOutput',0);

% basis transformation for Chern number
Qlistmat=cell2mat(Qlist');
Qshift1=mod(Qlistmat+qindex(1,:),1);
Qshift2=mod(Qlistmat+qindex(2,:),1);

Qshift1=mod(Qshift1,1);
Qshift2=mod(Qshift2,1);

perm1=zeros(size(Qshift1,1));
for i=1:length(Qshift1)
    j=find(abs(sum(abs(Qlistmat(i,:)-Qshift1).^2,2))<1e-5);
    perm1(i,j)=1;
end

perm2=zeros(size(Qshift2,1));
for i=1:length(Qshift2)
    j=find(abs(sum(abs(Qlistmat(i,:)-Qshift2).^2,2))<1e-5);
    perm2(i,j)=1;
end 

parameters.perm1=perm1;
parameters.perm2=perm2;

NQ=length(Qlist);
delta_tensor=zeros(NQ,NQ,NQ,NQ); %delta_tensor_{q_alpha,q_beta,q_gamma,q_delta}
for q_alpha_index=1:NQ
    for q_beta_index=1:NQ
        for q_gamma_index=1:NQ
            for q_delta_index=1:NQ
                qindex_alpha=Qlist{q_alpha_index};
                qindex_beta=Qlist{q_beta_index};
                qindex_gamma=Qlist{q_gamma_index};
                qindex_delta=Qlist{q_delta_index};
                deltafunc=qindex_gamma+qindex_delta-qindex_alpha-qindex_beta;
                delta_tensor(q_alpha_index,q_beta_index,q_delta_index,q_gamma_index)=all(abs(deltafunc-round(deltafunc))<1e-10);
            end
        end
    end
end
parameters.delta_tensor=delta_tensor;

end



