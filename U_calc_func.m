function U=U_calc_func(k,parameters)
% neighborlist{1}={[0,0]};
% neighborlist{2}={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]}; % direction on clock: 12, 10, 8, 6, 4, 2
% neighborlist{3}={[-1,-1],[1,-2],[2,-1],[1,1],[-1,2],[-2,1]}; %directin on clock: 11, 9, 7, 5, 4, 1
% neighborlist{4}={[-2,0],[0,-2],[2,-2],[2,0],[0,2],[-2,2]}; %direction on clock: 12, 10, 8, 6, 4, 2
% neighborlist{5}={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]}; %counterclockwise from the first point in Quadrant II
% neighborlist{6}={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]}; %direction on clock: 12, 10, 8, 6, 4, 2
% neighborlist{7}={[-4, 2], [-2, -2], [-2, 4], [2, -4], [2, 2], [4, -2]}; %direction TBD
% neighborlist{8}={[-4,1],[-4,3],[-3,-1],[-3,4],[-1,-3],[-1,4],[1,-4],[1,3],[3,-4],[3,1],[4,-3],[4,-1]};
% neighborlist{9}={[-4,0],[-4,4],[0,-4],[0,4],[4,-4],[4,0]};
% neighborlist{10}={[-5,2],[-5,3],[-3,-2],[-3,5],[-2,-3],[-2,5],[2,-5],[2,3],[3,-5],[3,2],[5,-3],[5,-2]};
% neighborlist{11}={[-5,1],[-5,4],[-4,-1],[-4,5],[-1,-4],[-1,5],[1,-5],[1,4],[4,-5],[4,1],[5,-4],[5,-1]};
% neighborlist{12}={[-5,0],[-5,5],[0,-5],[0,5],[5,-5],[5,0]};
% neighborlist{13}={[-6,3],[-3,-3],[-3,6],[3,-6],[3,3],[6,-3]};
% neighborlist{14}={[-6,2],[-6,4],[-4,-2],[-4,6],[-2,-4],[-2,6],[2,-6],[2,4],[4,-6],[4,2],[6,-4],[6,-2]};
% neighborlist{15}={[-6,1],[-6,5],[-5,-1],[-5,6],[-1,-5],[-1,6],[1,-6],[1,5],[5,-6],[5,1],[6,-5],[6,-1]};
% neighborlist{16}={[-6,0],[-6,6],[0,-6],[0,6],[6,-6],[6,0]};
% neighborlist{17}={[-7,3],[-7,4],[-4,-3],[-4,7],[-3,-4],[-3,7],[3,-7],[3,4],[4,-7],[4,3],[7,-4],[7,-3]};
% neighborlist{18}={[-7,3],[-7,4],[-4,-3],[-4,7],[-3,-4],[-3,7],[3,-7],[3,4],[4,-7],[4,3],[7,-4],[7,-3]};
% neighborlist{19}={[-7,1],[-7,6],[-6,-1],[-6,7],[-1,-6],[-1,7],[1,-7],[1,6],[6,-7],[6,1],[7,-6],[7,-1]};
% neighborlist{20}={[-7,0],[-7,7],[0,-7],[0,7],[7,-7],[7,0]};
% neighborlist{21}={[-8,3],[-8,5],[-7,0],[-7,7],[-5,-3],[-5,8],[-3,-5],[-3,8],[0,-7],[0,7],[3,-8],[3,5],[5,-8],[5,3],[7,-7],[7,0],[8,-5],[8,-3]};
% neighborlist{22}={[-8,2],[-8,6],[-6,-2],[-6,8],[-2,-6],[-2,8],[2,-8],[2,6],[6,-8],[6,2],[8,-6],[8,-2]};
% neighborlist{23}={[-8,1],[-8,7],[-7,-1],[-7,8],[-1,-7],[-1,8],[1,-8],[1,7],[7,-8],[7,1],[8,-7],[8,-1]};
% neighborlist{24}={[-8,0],[-8,8],[0,-8],[0,8],[8,-8],[8,0]};
neighborlist=generate_neighbor(40);

xrange=-5*sqrt(3)*parameters.aM:parameters.aM/20:5*sqrt(3)*parameters.aM;
yrange=-2*parameters.aM:parameters.aM/20:12*parameters.aM;
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

