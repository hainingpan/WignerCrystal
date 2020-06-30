function i=n2i(s,n)
%translate n position to linear index
i=intersect(find(n(:,1)==s(1)),find(n(:,2)==s(2)));
end