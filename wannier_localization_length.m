function loc=wannier_localization_length(parameters)
xrange=-3*parameters.aM:parameters.aM/20:3*parameters.aM;
yrange=-3*parameters.aM:parameters.aM/20:3*parameters.aM;
[rx,ry]=meshgrid(xrange,yrange);
[wbgrid,wtgrid]=w_rec({[0,0]},rx,ry,parameters);
psi2=abs(wbgrid(:,ceil(end/2))).^2;

% figure;
% surf(xrange/parameters.aM,yrange/parameters.aM,abs(wtgrid).^2,'edgecolor','none');
% view(2);
% axis tight;
% xlabel('x/a_M')
% ylabel('y/a_M')
% title('|\psi_t|^2')
% daspect([1,1,1])

yrange_q=-3*parameters.aM:parameters.aM/1000:3*parameters.aM;
wbabs2=interp1(yrange,psi2/max(psi2),yrange_q);
% figure;
% plot(yrange_q/parameters.aM,wbabs2);
% xlabel('y/a_M')
% ylabel('|\psi_b|^2 (normalized by the peak)')

[~,index]=min(abs(wbabs2-0.5));
loc=abs(yrange_q(index)/parameters.aM);
end