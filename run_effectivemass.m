thetalist=1:0.5:5;
clear nulist emlist
for i=1:length(thetalist)
    [nulist(:,:,i),emlist(:,:,i)]=effectivemass(thetalist(i));
end
kappalist=zeros(length(thetalist),1);
gammalist=zeros(length(thetalist),1);

for i=1:length(thetalist)
    kappa=emlist(nulist(:,:,i)==min(nulist(:,:,i),[],'all'),i);
    gamma=emlist(nulist(:,:,i)==max(nulist(:,:,i),[],'all'),i);
    kappalist(i)=1/sqrt(kappa);
    gammalist(i)=1/sqrt(gamma);
end

figure;
plot(thetalist,kappalist,'.-','displayname','\kappa');
hold on;
plot(thetalist,gammalist,'.-','displayname','\gamma');
legend;
xlabel('\theta');
ylabel('m^*/m_0');
