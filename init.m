function S=init(m,parameters)
%m=<..>=[SAx,SBx,SCx;SAy,SBy,SCy;SAz,SBz,SCz]
%Sr(a,b,c) a,b=1(up) or 2(down); c= 1(A), 2(B), 3(C)
Sr(1,1,:)=(1-m(3,:))/2;
Sr(1,2,:)=(m(1,:)-1i*m(2,:))/(-2);
Sr(2,1,:)=(m(1,:)+1i*m(2,:))/(-2);
Sr(2,2,:)=(m(3,:)+1)/2;
Q=[0,0;parameters.kn-parameters.kp;parameters.kp-parameters.kn];
r=[0,0;sqrt(3)/2*parameters.aM,1/2*parameters.aM;sqrt(3)*parameters.aM,parameters.aM];

S=zeros(2,2,3);
for i=1:3
    S(:,:,i)=1/3*((exp(-1i*Q(i,:)*r(1,:)'))*Sr(:,:,1)+...
        (exp(-1i*Q(i,:)*r(2,:)'))*Sr(:,:,2)+...
        (exp(-1i*Q(i,:)*r(3,:)'))*Sr(:,:,3));
end
end
