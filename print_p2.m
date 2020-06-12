param=mainTri2();
n=size(param.r,1);
p=zeros(n);
for i=1:n
    for j=i+1:n
        p(i,j)=pot2(param.uclist(i,:),param.uclist(j,:),param);
    end
end

p=p+p.';
p2=diag(2*sum(p))-p;
p2=int32(ceil(p2*1e7));
string="";
for i=1:n
    string=strcat(string,"%d ");
end
string=strcat(string,"\n");
file = fopen('p2.txt','w');
fprintf(file,string,p2);
fclose(file);
