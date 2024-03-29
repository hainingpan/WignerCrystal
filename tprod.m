function z=tprod(x,y,i,j,k,l)
%tensor product preserving certain index k,l, sum over i,j
%output index: first # of k follows the order of x, then goes the unsummed
%indices from x to y

xsize=size(x);
ysize=size(y);
if xsize(i)~=ysize(j)
    error('Must tensor product along the same dimensions');
end

if xsize(k)~=ysize(l)
    error('Must perserve along the same dimensions');
end

remdims_x = setdiff(1:ndims(x),[i,k]);
remdims_y = setdiff(1:ndims(y),[j,l]);
xorder=[i,k,remdims_x];
yorder=[j,l,remdims_y];
x=permute(x,xorder);
y=permute(y,yorder);
x=repmat(x,[xorder./xorder,ysize(remdims_y)]);
y=repmat(y,[yorder./yorder,xsize(remdims_x)]);
y=permute(y,[1:(length(j)+length(l)),...
    (length(j)+length(l)+length(remdims_y)+1):(length(j)+length(l)+length(remdims_y)+length(remdims_x)),...
    (length(j)+length(l))+1:(length(j)+length(l)+length(remdims_y))]);
if length(i)~=0
    z=squeeze(sum(x.*y,[1:length(i)]));
else
    z=squeeze(x.*y);
end

