function r_s=rs(epsilon,parameters)
n=parameters.nu(1)/parameters.nu(2)/((parameters.aM)^2*sqrt(3)/2);
r=1/sqrt(pi*n);

%4*pi*epsilon=e^2/(h*c*alpha)=1/alpha in natual unit
alpha=1/137;
rb=1/alpha*epsilon/(parameters.m);
r_s=r/rb;
end
