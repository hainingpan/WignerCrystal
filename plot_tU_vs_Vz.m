Ulist={};
% for i=1:length(Ustore)
%     Ulist{i}=arrayfun(@(x) Ustore{i}{x}(1),1:length(Ustore{i}));
% end
% for i=1:length(Ustore{1})
%     Ulist{i}=arrayfun(@(x) Ustore{x}{i}(1),1:length(Ustore));
% end
% figure;
% hold on;
% for i=1:length(Ulist)
%     plot(real(Ulist{i}));
% end
tlist={};
for i=1:length(tstore{1})
    tlist{i}=arrayfun(@(x) tstore{x}{i}(1),1:length(tstore));
end
figure;
hold on;
for i=2:4
    plot(Vzlist,cos(angle(tlist{i})));
end