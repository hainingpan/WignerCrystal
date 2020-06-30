function param=mainTri(varargin)
p=inputParser;
addParameter(p,'W',3);
addParameter(p,'L',3);
addParameter(p,'a',1);
% addParameter(p,'left',1);
% addParameter(p,'right',2);

parse(p,varargin{:});
a=p.Results.a;

param=struct('W',p.Results.W,'L',p.Results.L,'a',a);
param.a1=[sqrt(3)/2,1/2];
param.a2=[sqrt(3)/2,-1/2];
param.d=a*norm(param.a1);