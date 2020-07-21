function [sum,eigvec]=honeycomb_h2(phi,k,parameters)
a1=-parameters.aM1;
a2=parameters.aM2;
a3=-(a1+a2);
b1=a2-a3;
b2=a3-a1;
b3=a1-a2;

t2=.1;
ham=[2*t2*(cos(k*b1'+phi)+cos(k*b2'+phi)+cos(k*b3'+phi)),exp(1i*k*a1')+exp(1i*k*a2')+exp(1i*k*a3'),0;
    exp(-1i*k*a1')+exp(-1i*k*a2')+exp(-1i*k*a3'),2*t2*(cos(k*b1'-phi)+cos(k*b2'-phi)+cos(k*b3'-phi)),0;
    0,0,0];

ham=[0,1,0;0,0,1;1,0,0]'*ham*[0,1,0;0,0,1;1,0,0];
[eigvec,eigval]=eig(ham);
eigval=diag(eigval);
[~,I]=sort(eigval);
eigvec=eigvec(:,I);
sum=0; %alpha',alpha
for i=1:1
    sum=sum+conj(eigvec(:,i))*eigvec(:,i).';
end

end