function neighborlist=generate_neighbor(nshell)
a1=[0,-1];
a2=[sqrt(3)/2,-1/2];
counter=1;
for xindex=-nshell:nshell
    for yindex=max(-nshell,-nshell-xindex):min(nshell-xindex,nshell)
        x(counter)=xindex;
        y(counter)=yindex;
        r(counter)=norm(xindex*a1+yindex*a2)^2;
        counter=counter+1;
    end
end

[~,I]=sort(r);
xsorted=x(I);
ysorted=y(I);
brk=find(diff(r(I))>0.5);
first=[1,brk+1];
last=[brk,length(r)];
for i=1:length(first)
    coor=[xsorted(first(i):last(i));ysorted(first(i):last(i))]';
    clear neighbor
    for j=1:size(coor,1)
        neighbor{j}=coor(j,:);
    end        
    neighborlist{i}=neighbor;
end

    

        