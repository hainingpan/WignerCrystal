function f12=grad12(f,h1,h2)
f1=(-1/2*f([end,1:(end-1)],:)+1/2*f([2:end,1],:))/h1;
f12=(-1/2*f1(:,[end,1:(end-1)])+1/2*f1(:,[2:end,1]))/h2;

