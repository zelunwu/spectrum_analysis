function [R_tau,R_coef] = autocorr(ts,dt,unit)
% [R_tau,R_coef] = AUTOCORR(ts,dt,unit)
% A function to calculate Autocovariance
% Input:
% 	ts: time series
% 	dt: sampling rate
%
% Output:
% 	R_tau: Auto-covariance
%	R_coef: Auto-correlation, R_coef = R_tau/var(ts)
% 
%-------------------------------------------------------------------------%
%   Author: Zelun Wu,
%   Ph.D. student of Physical Oceanography
%   University of Delaware, Xiamen University
%   zelunwu@udel.edu, zelunwu@stu.xmu.edu.cn
%-------------------------------------------------------------------------%


N = length(ts);
T = (N-1) * dt; % total time

for k = 0:N-1
    sum_k = 0;
    for n = 1:N-k
        sum_k = sum_k+ ts(n)*ts(n+k);
    end
    R_tau(k+1) = sum_k;
end
R_tau = R_tau./N;
R_coef = R_tau/var(ts);

%% Plot
time = reshape([0:N-1]*dt,size(ts)); %time
plot(time,R_coef,'linewidth',2);
xlabel(['\tau (lag, ',unit,')']);
ylabel('Autocorrelation Coefficient');
set(gca,'fontsize',20);

end