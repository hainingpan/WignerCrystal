function re=hubbardU_fft(wbgrid1,wtgrid1,wbgrid2,wtgrid2,rx,ry,parameters)
%Obsolete, hubbardU_fft_2 instead; for record purpose only
[Nx,Ny]=size(rx);
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
w2=reshape(wb1k.*wb2k+wt1k.*wt2k,Nx,Ny);

F=griddedInterpolant(kxmap',kymap',w2);
func=@(x,y) alpha/(2*pi)*1./(sqrt(x.^2+y.^2)).*(F(x',y')').*(1-exp(-parameters.d*(sqrt(x.^2+y.^2))));
re=quad2d(func,kxlist(1),kxlist(end),kylist(1),kylist(end),'FailurePlot',false,'MaxFunEvals',5000);

end