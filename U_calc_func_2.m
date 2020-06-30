function U=U_calc_func_2(k,parameters)
%use shifted wannier state
neighborlist=generate_neighbor(40);

xrange=-3*parameters.aM:parameters.aM/20:3*parameters.aM;
yrange=-3*parameters.aM:parameters.aM/20:3*parameters.aM;
[rx,ry]=meshgrid(xrange,yrange);


% [wbgrid,wtgrid]=w_rec([neighborlist{1:k+1}],rx,ry,parameters);
neighborlist2=cellfun(@(x)x{1},neighborlist,'UniformOutput',false);
[wbgrid,wtgrid]=w_rec([neighborlist2(1:1)],rx,ry,parameters);


U={};
Uint=hubbardU_fft_2(wbgrid,wtgrid,[neighborlist2(1:k+1)],rx,ry,parameters);
for i=1:k+1
%     Uint=hubbardU_fft(wbgrid,wt{1},wb{i},wt{i},rx,ry,parameters);
    for j=1:length(neighborlist{i})
        U{i}(j)=Uint(i);
    end
end

end

