function en=pot(x,y,d,param)

r=dist(x,y,param);
en=1/r-1/sqrt(r^2+d^2);
end