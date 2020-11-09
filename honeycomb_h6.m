function [sum,eigvec2,eigval]=honeycomb_h6(phi,k,parameters)
a1=-parameters.aM1;
a2=parameters.aM2;
a3=-(a1+a2);
b1=a2-a3;
b2=a3-a1;
b3=a1-a2;

t2=.1;
ham11=[2*t2*(cos(k*b1'+phi)+cos(k*b2'+phi)+cos(k*b3'+phi)),exp(1i*k*a1')+exp(1i*k*a2')+exp(1i*k*a3'),0;
    exp(-1i*k*a1')+exp(-1i*k*a2')+exp(-1i*k*a3'),2*t2*(cos(k*b1'-phi)+cos(k*b2'-phi)+cos(k*b3'-phi)),0;
    0,0,10];
ham22=10*eye(3);
ham22(3,3)=-10;
ham=[ham11,zeros(3);zeros(3),ham22];
[eigvec,eigval]=eig(ham);
eigval=diag(eigval);
eigvec2=eigvec;
[~,I]=sort(eigval);
eigvec=eigvec(:,I);
sum=0; %alpha',alpha
for i=1:2
    sum=sum+conj(eigvec(:,i))*eigvec(:,i).';
end

end