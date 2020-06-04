function [val,vec]=energyTMD(kx,ky,parameters)
%energy at (kx,ky)
DeltaTmat=parameters.DeltaTmat;
DeltaTTmat=parameters.DeltaTTmat;
Deltabmat=parameters.Deltabmat;
Deltatmat=parameters.Deltatmat;
m=parameters.m;
kp=parameters.kp;
kn=parameters.kn;
bM1=parameters.bM1;
bM2=parameters.bM2;
Vz=parameters.Vz;
h1index=parameters.h1index;
h2index=parameters.h2index;
kplist=[kx,ky]-kp+h1index(:)*bM1+h2index(:)*bM2;
knlist=[kx,ky]-kn+h1index(:)*bM1+h2index(:)*bM2;
H11=-1/(2*m)*diag(dot(kplist,kplist,2))+Deltabmat+Vz/2*eye(length(Deltabmat));
H22=-1/(2*m)*diag(dot(knlist,knlist,2))+Deltatmat-Vz/2*eye(length(Deltabmat));
H12=DeltaTmat;
H21=DeltaTTmat;
H=[H11,H12;H21,H22];
[vec,val]=eig(H);
val=diag(val);
val=val(end:-1:1);
vec=vec(:,end:-1:1);
end


