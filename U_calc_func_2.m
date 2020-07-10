function U=U_calc_func_2(k,parameters)
%use shifted wannier state
neighborlist=generate_neighbor(101);

xrange=-3*parameters.aM:parameters.aM/20:3*parameters.aM;
yrange=-3*parameters.aM:parameters.aM/20:3*parameters.aM;
[rx,ry]=meshgrid(xrange,yrange);


% [wbgrid,wtgrid]=w_rec([neighborlist{1:k+1}],rx,ry,parameters);
neighborlist2=cellfun(@(x)x{1},neighborlist,'UniformOutput',false);
[wbgrid,wtgrid]=w_rec(neighborlist2(1),rx,ry,parameters);
neighbordist=cellfun(@(x)norm(x*[parameters.aM1;parameters.aM2]),neighborlist2);
k1=find(neighbordist>2*parameters.d,1);
if isempty(k1)
    k1=k+1;
end

U={};
Uint=hubbardU_fft_2(wbgrid,wtgrid,[neighborlist2(1:k1)],rx,ry,parameters);
alpha=0.00729735;
Uint2=alpha./(neighbordist(k1+1:k+1))-alpha./sqrt(neighbordist(k1+1:k+1).^2+(parameters.d)^2);
Utot=[Uint,Uint2];
for i=1:k+1
%     Uint=hubbardU_fft(wbgrid,wt{1},wb{i},wt{i},rx,ry,parameters);
    for j=1:length(neighborlist{i})
        U{i}(j)=Utot(i);
    end
end

end

