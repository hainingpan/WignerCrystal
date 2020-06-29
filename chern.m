function berrycur=chern(wfall,parameters)
n=sqrt(size(wfall,1));
NL=size(wfall,2);
% bM1=parameters.bM1;
% bM2=parameters.bM2;
% kp=parameters.kp;
% kn=parameters.kn;
% a1=bM2/(2*n);
% a2=(-2*bM1-bM2)/2/(2*n);
% omega=abs(cross([bM1,0],[bM2,0]));
% omega=omega(3);

level=1:NL*parameters.nu(1)/(2*parameters.nu(2));

% xrange=-n:n;
% yrange=-n:n;
% kxmap=zeros(2*n+1,2*n+1);
% kymap=zeros(2*n+1,2*n+1);
% kcx2map=zeros(2*n,2*n); %after shifting to Hexagon
% kcy2map=zeros(2*n,2*n);  %after shifting to Hexagon
% umap=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
% omega=abs(cross([a1,0],[a2,0]));
% omega=omega(3);


% Nx=length(xrange);
% Ny=length(yrange);
% for xindex=1:Nx
%     kx=xrange(xindex);
%     for yindex=1:Ny
%         ky=yrange(yindex);
%         k=kx*a1+ky*a2;
% %         [~,vec]=energyTMD(k(1),k(2),parameters);
%         umap(xindex,yindex,:)=vec(:,level);
%         kxmap(xindex,yindex)=k(1);
%         kymap(xindex,yindex)=k(2);
%     end
% end
% kxmap=reshape(kxlist,[n,n]);
% kymap=reshape(kylist,[n,n]);

wfmap=reshape(wfall,[n,n,size(wfall,2),size(wfall,3)]); %Nkx,Nky,level index, wf component
% kcxmap=(kxmap(1:end-1,1:end-1)+kxmap(1:end-1,2:end)+kxmap(2:end,1:end-1)+kxmap(2:end,2:end))/4;
% kcymap=(kymap(1:end-1,1:end-1)+kymap(1:end-1,2:end)+kymap(2:end,1:end-1)+kymap(2:end,2:end))/4;



% for xindex=1:Nx-1
%     for yindex=1:Ny-1
%         k=[kcxmap(xindex,yindex),kcymap(xindex,yindex)];
%         shift=[0,0];
%         if (k(1)<=0) && (k(2)>=2*kn(2)/bM1(1)*k(1)+2*kn(2))
%             shift=a1*2*n;
%         end
%         if (k(1)>=0) && (k(2)>=-2*kn(2)/bM1(1)*k(1)+2*kn(2))
%             shift=-a2*2*n;
%         end
%         if (k(1)<=0) && (k(2)<=-2*kn(2)/bM1(1)*k(1)-2*kn(2))
%             shift=a2*2*n;
%         end
%         if (k(1)>=0) && (k(2)<=2*kn(2)/bM1(1)*k(1)-2*kn(2))
%             shift=-a1*2*n;
%         end   
%         kcx2map(xindex,yindex)=k(1)+shift(1);
%         kcy2map(xindex,yindex)=k(2)+shift(2);
%     end
% end
    
    
    % kcx2map=(kx2map(1:end-1,1:end-1)+kx2map(1:end-1,2:end)+kx2map(2:end,1:end-1)+kx2map(2:end,2:end))/4;
% kcy2map=(ky2map(1:end-1,1:end-1)+ky2map(1:end-1,2:end)+ky2map(2:end,1:end-1)+ky2map(2:end,2:end))/4;

% for i=1:NL
%     umap=squeeze(wfmap(:,:,i,:));
%     bcmap(:,:,i)=-angle(dot((umap(1:end-1,1:end-1,:)),umap(1:end-1,2:end,:),3)...
%     .*dot((umap(1:end-1,2:end,:)),umap(2:end,2:end,:),3)...
%     .*dot((umap(2:end,2:end,:)),umap(2:end,1:end-1,:),3)...
%     .*dot((umap(2:end,1:end-1,:)),umap(1:end-1,1:end-1,:),3));    
% end

    umap=(wfmap(:,:,level,:));
%     u12=ttt2(conj(umap(1:end-1,1:end-1,:,:)),umap(1:end-1,2:end,:,:),[4],[4],[1,2],[1,2]);
%     u23=ttt2(conj(umap(1:end-1,2:end,:,:)),umap(2:end,2:end,:,:),[4],[4],[1,2],[1,2]);
%     u34=ttt2(conj(umap(2:end,2:end,:,:)),umap(2:end,1:end-1,:,:),[4],[4],[1,2],[1,2]);
%     u41=ttt2(conj(umap(2:end,1:end-1,:,:)),umap(1:end-1,1:end-1,:,:),[4],[4],[1,2],[1,2]);
%     prod123=ttt2(u12,u23,[4],[3],[1,2],[1,2]);
%     prod1234=ttt2(prod123,u34,[4],[3],[1,2],[1,2]);
%     prod12341=ttt2(prod1234,u41,[4],[3],[1,2],[1,2]);
% [N1,N2,~,~]=size(prod12341);
berrycur=zeros(n-1);
for i=1:n
    for j=1:n
        u1=squeeze(umap(i,j,:,:));
        
        u2=squeeze(umap(mod(i,n)+1,j,:,:))*((i+1>n)*kron(eye(4),[0,1;1,0])+(i+1<=n)*eye(8));
        
        u3=squeeze(umap(mod(i,n)+1,mod(j,n)+1,:,:))*((i+1>n)*kron(eye(4),[0,1;1,0])+(i+1<=n)*eye(8))*...
            ((j+1>n)*kron(eye(2),kron([0,1;1,0],eye(2)))+(j+1<=n)*eye(8));
        
        u4=squeeze(umap(i,mod(j,n)+1,:,:))*((j+1>n)*kron(eye(2),kron([0,1;1,0],eye(2)))+(j+1<=n)*eye(8));

        U12=conj(u1)*u2.';
        U23=conj(u2)*u3.';
        U34=conj(u3)*u4.';
        U41=conj(u4)*u1.';
        berrycur(i,j)=angle(det(U12*U23*U34*U41));
    end
end

% eiv=zeros(n,2);
% for i=1:n
%     prod=eye(length(level));
%     for j=1:n
%         
%         if (j+1)~=n+1
%             prod=prod*squeeze(conj(umap(i,j,:,:)))*squeeze(umap(i,j+1,:,:)).';
%         else
%             check=squeeze(conj(umap(i,j,:,:)))*kron(eye(4),[0,1;1,0])*squeeze(umap(i,1,:,:)).';
%             disp(check);
%             prod=prod*check;
%         end
%         eigval=eig(prod);
%         eiv(i,:)=angle(eigval)/(pi);
%     end
% end

% ch=sum(berrycur(:))/(2*pi);
% ch=squeeze(sum(bcmap,[1,2]))/(2*pi);
end