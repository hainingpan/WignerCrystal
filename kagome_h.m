function [sum,eigvec]=kagome_h(phi,s1,k,parameters)
a1=-parameters.aM1;
a2=parameters.aM2;
a3=-(a1+a2);
% ham=2*[0,exp(1i*phi)*cos(k*a1'),exp(-1i*phi)*cos(k*a3'),0;
%        exp(-1i*phi)*cos(k*a1'),0,exp(1i*phi)*cos(k*a2'),0;
%        exp(1i*phi)*cos(k*a3'),exp(-1i*phi)*cos(k*a2'),0,0;
%        0,0,0,10];
ham=2*[0,exp(1i*phi)*(exp(-1i*k*a1')+s1*exp(1i*k*a1')),exp(-1i*phi)*(exp(1i*k*a3')+s1*exp(-1i*k*a3')),0;
       exp(-1i*phi)*(exp(-1i*k*a1')+s1*exp(1i*k*a1')),0,exp(1i*phi)*(exp(-1i*k*a2')+s1*exp(1i*k*a2')),0;
       exp(1i*phi)*(exp(-1i*k*a3')+s1*exp(1i*k*a3')),exp(-1i*phi)*(exp(1i*k*a2')+s1*exp(-1i*k*a2')),0,0;
       0,0,0,10];
[eigvec,eigval]=eig(ham);
eigval=diag(eigval);
[~,I]=sort(eigval);
eigvec=eigvec(:,I);
sum=0; %alpha',alpha
for i=1:2
    sum=sum+conj(eigvec(:,i))*eigvec(:,i).';
end

end