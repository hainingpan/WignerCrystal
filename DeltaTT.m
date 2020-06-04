function re=DeltaTT(h1,h2,parameters)
%Delta_T^dagger in Fourier space
w=parameters.w;
re=w*((h1==0).*(h2==0)+(h1==1).*(h2==1)+(h1==1).*(h2==0));
end