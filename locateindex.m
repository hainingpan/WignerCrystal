function i=locateindex(Q,Qlistmod)
Q=[mod(Q(1),1),mod(Q(2),1)];
for i=1:length(Qlistmod)
    if isequal(Q,Qlistmod{i})
        return;
    end
end
i=-1;
