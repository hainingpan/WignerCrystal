hold on;
x=generate_neighbor(4);
for xindex=1:length(x)
    for yindex=1:length(x{xindex})
        pts=(x{xindex}{yindex});
        cor=pts(2)*[0,-1]+pts(1)*[sqrt(3)/2,-1/2];
        scatter3(cor(1),cor(2),1e4,100,'r.')
    end
end
axis([-3,3,-3,3])