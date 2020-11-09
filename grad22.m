function f22=grad22(f,h)
f22=(f(:,[2:end,1])+f(:,[end,1:(end-1)])-2*f)/h^2;
