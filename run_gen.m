% for i=3:0.04:3.96
%     disp(i);
%     sweep_ep_theta_gen([1,2],[0],[i]);
% end

sweep_ep_theta([1,6],1:60,[1],0,'phase1,6_theta(1.00,1.00,1)_Vz(0.0).mat');
% sweep_ep_theta([2,8],1:60,[1],0,'phase1,4_theta(1.00,1.00,1)_Vz(0.0).mat');

% 
% sweep_ep_theta([5,12]*2,1:35,[4],0,'phase10,24_theta(4.00,4.00,1)_Vz(0.0).mat');
% sweep_ep_theta([5,12]*5,1:35,[4],0,'phase10,24_theta(4.00,4.00,1)_Vz(0.0).mat');
% sweep_ep_theta([5,12]*6,1:35,[4],0,'phase10,24_theta(4.00,4.00,1)_Vz(0.0).mat');
