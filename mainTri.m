function param=mainTri(varargin)
p=inputParser;
addParameter(p,'W',3);
addParameter(p,'L',3);
addParameter(p,'phi',1);
% addParameter(p,'left',1);
% addParameter(p,'right',2);

parse(p,varargin{:});
phi=p.Results.phi;
phi=[phi,-phi,phi,-phi,phi,-phi];
Jx=cos(2*phi);
Jy=cos(2*phi);
Jz=ones(1,6);
Z=sin(2*phi);
% param=struct('W',p.Results.W,'L',p.Results.L,'phi',phi,'left',p.Results.left,'right',p.Results.right,'Jx',Jx,'Jy',Jy,'Jz',Jz,'Z',Z);
param=struct('W',p.Results.W,'L',p.Results.L,'phi',phi,'Jx',Jx,'Jy',Jy,'Jz',Jz,'Z',Z);
param.a1=[sqrt(3)/2,1/2];
param.a2=[sqrt(3)/2,-1/2];