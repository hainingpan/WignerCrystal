function re=screening(x,y,parameters)
if parameters.s==inf
    re=1;
else
    re=(parameters.s*sqrt(x.^2+y.^2))./(sqrt(1+parameters.s^2*(x.^2+y.^2)));
end