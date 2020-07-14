for i=1:size(energyall,2)
    wftmp=abs(squeeze(wfall(1,i,:))).^2;
    [~,idx]=max([sum(wftmp(1:end/2)),sum(wftmp(end/2+1:end))]);
    spinpolorize(i)=idx;
end