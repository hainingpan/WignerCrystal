function re=tb(bond,t,kxmap,kymap,parameters)
assert(length(bond)==length(t),'Length mismatch');
n=length(bond);
aM1=parameters.aM1;
aM2=parameters.aM2;
re=0;
for i=1:n
    a=bond{i}(1)*aM1+bond{i}(2)*aM2;
    re=re+exp(-1i*(kxmap*a(1)+kymap*a(2)))*t(i);
end