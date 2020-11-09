function f11=grad11(f,h)
f11=(f([2:end,1],:)+f([end,1:(end-1)],:)-2*f)/h^2;
