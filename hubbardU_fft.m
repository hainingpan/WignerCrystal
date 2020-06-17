function re=hubbardU_fft(wbgrid1,wtgrid1,wbgrid2,wtgrid2,rx,ry,parameters)
% bM1=parameters.bM1;
% bM2=parameters.bM2;
% [rx,ry]=meshgrid(linspace(-2*sqrt(3)*parameters.aM,2*sqrt(3)*parameters.aM,101));
[Nx,Ny]=size(rx);
% [wbgrid,wtgrid]=w(R,rx,ry,parameters);
alpha=0.00729735; %e^2/(4*pi*epsilon0); c.f. Quicks Notes on oneNotes

Lx=rx(1,end)-rx(1,1);
Ly=ry(end,1)-ry(1,1);

kxlist=2*pi/Lx*(-floor(Nx/2):floor((Nx-1)/2));
kylist=2*pi/Ly*(-floor(Ny/2):floor((Ny-1)/2));

wb1k=fftshift(fft2(abs(wbgrid1).^2))*Lx*Ly/(Nx*Ny);
wt1k=fftshift(fft2(abs(wtgrid1).^2))*Lx*Ly/(Nx*Ny);

wb2k=fftshift(fft2(abs(wbgrid2).^2))*Lx*Ly/(Nx*Ny);
wb2k=reshape(wb2k(end:-1:1),Nx,Ny);
wt2k=fftshift(fft2(abs(wtgrid2).^2))*Lx*Ly/(Nx*Ny);
wt2k=reshape(wt2k(end:-1:1),Nx,Ny);

[kxmap,kymap]=meshgrid(kxlist,kylist);
% Nkx=200;
% Nky=200;
% kxlist=linspace(-3*bM2(1),3*bM2(1),Nkx);
% kylist=linspace(-3*bM2(1),3*bM2(1),Nky);

% [kxmap,kymap]=meshgrid(kxlist,kylist);
% Fb=exp(1i*(repmat(reshape(kxmap,1,1,[]),Nx,Ny,1).*rx+repmat(reshape(kymap,1,1,[]),Nx,Ny,1).*ry)).*abs(wbgrid).^2;
% intb=trapz(rx(1,:),trapz(ry(:,1),Fb,2));
% Ft=exp(1i*(repmat(reshape(kxmap,1,1,[]),Nx,Ny,1).*rx+repmat(reshape(kymap,1,1,[]),Nx,Ny,1).*ry)).*abs(wtgrid).^2;
% intt=trapz(rx(1,:),trapz(ry(:,1),Ft,2));
% w2=reshape(abs(wbk).^2+abs(wtk).^2,Nx,Ny);
w2=reshape(wb1k.*wb2k+wt1k.*wt2k,Nx,Ny);

F=griddedInterpolant(kxmap',kymap',w2);
func=@(x,y) alpha/(2*pi)*1./(sqrt(x.^2+y.^2)).*(F(x',y')');
re=quad2d(func,kxlist(1),kxlist(end),kylist(1),kylist(end),'FailurePlot',false,'MaxFunEvals',5000);
% Fk=alpha/(2*pi)*1./(sqrt(kxmap.^2+kymap.^2)).*reshape(abs(wbk).^2+abs(wtk).^2,Nx,Ny);
% re=trapz(kxlist,trapz(kylist,Fk,2));
end