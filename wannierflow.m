function eiv=wannierflow(wfall,parameters)
%Wannier center flow; c.f. yu2011equivalent
n=sqrt(size(wfall,1));
NL=size(wfall,2);
level=1:NL*parameters.nu(1)/(2*parameters.nu(2));
wfmap=reshape(wfall,[n,n,size(wfall,2),size(wfall,3)]); %Nkx,Nky,level index, wf component
umap=(wfmap(:,:,level,:));

eiv=zeros(n,2);
for i=1:n
    prod=eye(length(level));
    for j=1:n
        
        if (j+1)~=n+1
            prod=prod*reshape(conj(umap(i,j,:,:)),[level(end),NL])*reshape(umap(i,j+1,:,:),[level(end),NL]).';
        else
%             check=squeeze(conj(umap(i,j,:,:)))*kron(eye(4),[0,1;1,0])*squeeze(umap(i,1,:,:)).';
            check=reshape(conj(umap(i,j,:,:)),[level(end),NL])*kron(eye(2),parameters.perm1)*reshape(umap(i,1,:,:),[level(end),NL]).';
            prod=prod*check;
        end
        eigval=eig(prod);
        eiv(i,:)=angle(eigval)/(pi);
    end
end
figure;plot(eiv,'.');
ylabel('\phi/\pi');
xlim([1,sqrt(size(wfall,1))]);
end