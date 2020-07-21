wfstore=zeros(size(wfall,3));
for i=1:size(energyall,2)
    wftmp=abs(squeeze(wfall(43,i,:))).^2;
    wfstore(:,i)=wftmp;
    [~,idx]=max([sum(wftmp(1:end/2)),sum(wftmp(end/2+1:end))]);
    spinpolorize(i)=idx;
end