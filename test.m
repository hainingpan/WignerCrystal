% ztmp=zeros(2,4,5,6,7);
% for i=1:2
%     for j=1:4
%         ztmp(i,j,:,:,:)=ttt(tensor(squeeze(x(i,j,:,:,:))),tensor(squeeze(y(:,:,j,i))),[1],[1]);
%     end
% end

% ave=average_init(parameters);

% [vec,val]=eig(h2);
%     val=real(diag(val));
%     [val,I]=sort(val);
%     vec=vec(:,I);
%     energyall2(k_beta_index,:)=val;
%     for ii=1:2*NQ
%         wfall2(k_beta_index,ii,:)=vec(:,ii);
%     end
% for i=1:36
% surf(xrange,yrange,abs(wb{i}),'edgecolor','none');view(2);title(num2str(norm(neighborlist2{i}*[0,-1;[sqrt(3)/2,-1/2]])^2))
% end
% for i=11:20
%     disp(i);
% tic;sweep_nu([1,3],10,10,length(generate_neighbor(i+4)));tt(i)=toc;
% end
sweep_nu([1,1],10,10);
sweep_nu([1,2],10,10);
sweep_nu([1,3],10,10);
sweep_nu([2,3],10,10);
sweep_nu([1,4],10,10);
sweep_nu([3,4],10,10);


