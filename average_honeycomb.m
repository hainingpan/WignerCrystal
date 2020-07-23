function [ave3,V2ave]=average_honeycomb(kxlist,kylist,parameters)
%calculate <c_{k+q_1,up}^dagger c_{k+q_2,up}>

Nk=length(kxlist);
ave=zeros(3,3,Nk); %alpha',alpha,k
wf=zeros(Nk,3,3);
for kindex=1:Nk
    kx=kxlist(kindex);
    ky=kylist(kindex);
    [ave(:,:,kindex),wf(kindex,:,:)]=honeycomb_h(pi/2,[kx,ky],parameters);
end

r={[0,0];[-parameters.aM1];[parameters.aM2-parameters.aM1]};
Q=parameters.Q;
NQ=length(Q);
expalphaQ=zeros(3,NQ); %alpha',q1
for alphaindex=1:3
    for Qindex=1:NQ
        expalphaQ(alphaindex,Qindex)=exp(1i*r{alphaindex}*Q{Qindex}');
    end
end

expalphaQ2=conj(expalphaQ);

prod1=ttt(tensor(permute(ave,[3,1,2])),tensor(expalphaQ2),[3],[1]);
prod2=ttt(tensor(permute(prod1,[1,3,2])),tensor(expalphaQ),[3],[1]);

ave2=permute(prod2,[1,3,2])/3; %k,q1,q2

ave3=zeros([size(ave2),2,2]); 
ave3(:,:,:,1,1)=ave2;
if parameters.nu==[4,6]
%     for kindex=1:Nk
%         kx=kxlist(kindex);
%         ky=kylist(kindex);
%         [ave(:,:,kindex),wf(kindex,:,:)]=honeycomb_h2(-pi/2,[kx,ky],parameters);
%     end
%     prod1=ttt(tensor(permute(ave,[3,1,2])),tensor(expalphaQ2),[3],[1]);
%     prod2=ttt(tensor(permute(prod1,[1,3,2])),tensor(expalphaQ),[3],[1]);
% 
%     ave2=permute(prod2,[1,3,2])/3; %k,q1,q2
%     ave3(:,:,:,2,2)=ave2;
end

V2=parameters.V2;

V2ave=ttt2(tensor(V2),tensor(ave3),[1],[1],[3],[2]); %q_alpha,k_beta,q_delta,q_gamma,sigma,sigma'

end
        

