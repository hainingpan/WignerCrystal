function re=Deltal(h1,h2,l,parameters)
%l=1 for bottom layer
%l=-1 for top layer
V=parameters.V;
psi=parameters.psi;
re=V*((h1==0).*(h2==-1)*exp(1i*l*psi)+...    %G1
    (h1==0).*(h2==1)*exp(-1i*l*psi)+...    %G4
    (h1==-1).*(h2==0)*exp(1i*l*psi)+...  `   %G5
    (h1==1).*(h2==0)*exp(-1i*l*psi)+...    %G2
    (h1==1).*(h2==1)*exp(1i*l*psi)+...    %G3
    (h1==-1).*(h2==-1)*exp(-1i*l*psi));       %G6
end