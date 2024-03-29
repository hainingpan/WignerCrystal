function sum=pot_symb(x,y,param)
sum=0;
nn=length(param.neighbor)-1;
U=sym('U',[1,nn]);
for i=1:size(param.Aneighbor2,1)
    y2=y+param.Aneighbor2(i,:);
    r=norm((y2-x)*[param.a1;param.a2])^2;
    index=find(abs(r-param.neighbor(2:end))<1e-2);
    if ~isempty(index) & index<=nn
        sum=U(index)+sum;
    end
end
end