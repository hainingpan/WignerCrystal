param=mainTri('W',3,'L',6);
% param=mainTri('W',2,'L',ceil(2));
param=construct_triangular(param);
n=length(param.r);
p=zeros(n);
a=0.1*norm(param.a1);
for i=1:n
    for j=i+1:n
        p(i,j)=pot(param.r(i,:),param.r(j,:),a,param);
    end
end

p=p+p.';
p2=diag(2*sum(p))-p;
p2=int32(ceil(p2*3e5));
string="";
for i=1:n
    string=strcat(string,"%d ");
end
string=strcat(string,"\n");
file = fopen('p2.txt','w');
fprintf(file,string,p2);
fclose(file);
