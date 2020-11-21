%find the common multiplier for n
function n=cm(delta)
delta=round(delta); 
nlist=3*find(arrayfun(@(k)mod(3^3*k^2,delta),1:20)==0);
n1=nlist(nlist>=21);
n=n1(1);
end