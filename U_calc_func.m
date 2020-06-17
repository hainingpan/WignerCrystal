function U=U_calc_func(k,parameters)
neighborlist{1}={[0,0]};
neighborlist{2}={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]}; % direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{3}={[-1,-1],[1,-2],[2,-1],[1,1],[-1,2],[-2,1]}; %directin on clock: 11, 9, 7, 5, 4, 1
neighborlist{4}={[-2,0],[0,-2],[2,-2],[2,0],[0,2],[-2,2]}; %direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{5}={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]}; %counterclockwise from the first point in Quadrant II
neighborlist{6}={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]}; %direction on clock: 12, 10, 8, 6, 4, 2
neighborlist{7}={[-4, 2], [-2, -2], [-2, 4], [2, -4], [2, 2], [4, -2]}; %direction TBD
neighborlist{8}={[-4,1],[-4,3],[-3,-1],[-3,4],[-1,-3],[-1,4],[1,-4],[1,3],[3,-4],[3,1],[4,-3],[4,-1]};
neighborlist{9}={[-4,0],[-4,4],[0,-4],[0,4],[4,-4],[4,0]};

xrange=-2*sqrt(3)*parameters.aM:parameters.aM/20:2*sqrt(3)*parameters.aM;
yrange=-2*parameters.aM:parameters.aM/20:6*parameters.aM;
[rx,ry]=meshgrid(xrange,yrange);


% [wbgrid,wtgrid]=w_rec([neighborlist{1:k+1}],rx,ry,parameters);
neighborlist2=cellfun(@(x)x{1},neighborlist,'UniformOutput',false);
[wbgrid,wtgrid]=w_rec([neighborlist2(1:k+1)],rx,ry,parameters);

counter=1;
for i=1:k+1
        wb{i}=wbgrid(:,:,counter);
        wt{i}=wtgrid(:,:,counter);
    counter=counter+1;
end
U={};
for i=1:k+1
    Uint=hubbardU_fft(wb{1},wt{1},wb{i},wt{i},rx,ry,parameters);
    for j=1:length(neighborlist{i})
        U{i}(j)=Uint;
    end
end

end

