function param=mainTri2()

param.d=0.1;
param.a1=[0,-1];
param.a2=[sqrt(3)/2,-1/2];

param.neighborlist{1}={[0,0]};
param.neighborlist{2}={[-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]}; % direction on clock: 12, 10, 8, 6, 4, 2
param.neighborlist{3}={[-1,-1],[1,-2],[2,-1],[1,1],[-1,2],[-2,1]}; %directin on clock: 11, 9, 7, 5, 4, 1
param.neighborlist{4}={[-2,0],[0,-2],[2,-2],[2,0],[0,2],[-2,2]}; %direction on clock: 12, 10, 8, 6, 4, 2
param.neighborlist{5}={[-2,-1],[-1,-2],[1,-3],[2,-3],[3,-2],[3,-1],[2,1],[1,2],[-1,3],[-2,3],[-3,2],[-3,1]}; %counterclockwise from the first point in Quadrant II
param.neighborlist{6}={[-3,0],[0,-3],[3,-3],[3,0],[0,3],[-3,3]}; %direction on clock: 12, 10, 8, 6, 4, 2

n=10;
counter=1;
N=3*n^2+3*n+1;
xlist=zeros(N,1);
ylist=zeros(N,1);
rxlist=zeros(N,1);
rylist=zeros(N,1);
for xindex=-n:n
    for yindex=max(-n,-n-xindex):min(n-xindex,n)
        xlist(counter)=xindex;
        ylist(counter)=yindex;
        r=xindex*param.a1+yindex*param.a2;
        rxlist(counter)=r(1);
        rylist(counter)=r(2);
        counter=counter+1;
    end
end

param.xlist=xlist;
param.ylist=ylist;
param.rxlist=rxlist;
param.rylist=rylist;

distsquare=sort(rxlist.^2+rylist.^2);
param.neighbor=distsquare(logical([1;diff(distsquare)>.5]));


% for i=1:length(param.neighborlist)
%     param.neighbor(i)=norm(param.neighborlist{i}{1}(1)*param.a1+param.neighborlist{i}{1}(2)*param.a2);    
% end
% param.Ulist=1./param.neighbor(2:end);

%% hexgonal supercell with 12 sites
param.A1=[2,2];
param.A2=[-2,4];
param.uclist=[[0,0];[-1,0];[0,-1];[1,-1];[1,0];[0,1];[-1,1];...
    [1,1];[0,2];[-1,2];[-2,2];[-2,1]];
param.r=param.uclist*[param.a1;param.a2];
% Aneighbor=[[0,0];[1,0];[0,1];[-1,1];[-1,0];[0,-1];[1,-1]];
Aneighbor=reshape([param.neighborlist{1}{:},param.neighborlist{2}{:},param.neighborlist{3}{:}],[2,13])';
%% rectangular supercell with 4 sites
% param.A1=[-1,2];
% param.A2=[-2,0];
% param.uclist=[[0,0];[-1,1];[-1,0];[-2,1]];
% param.r=param.uclist*[param.a1;param.a2];
% [X,Y]=meshgrid(-2:2,-2:2);
% Aneighbor=[X(:),Y(:)];


Aneighbor2=Aneighbor*[param.A1;param.A2];
param.Aneighbor2=Aneighbor2;


end