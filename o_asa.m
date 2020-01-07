function [apsd,freq] = o_asa(ts,dt,M_seg,option)
%-------------------------------------------------------------------------%
% Function:
% [apsd,freq] = o_asa(ts,dt,M_seg,option)
% Auto Spectrum Analysis
% Based on the Discrete Fourier Transform function: o_dft.m.
% 
% f(n+1) = n*dt/N_segment, n = 0,1,2,...,N_segment-1
% apsd = 1/Time_segment*|Xf|.^2, where Xf is the Fourier Transform of ts.
% If M_seg>1, then apsd is the mean apsd of M_segments divided from ts.
% 
% Input:
%     ts: Input time series. shoule be an N*1 or 1*N sequence.
%     dt: Sampling rate in time domain. E.g., 1 hour.
%     M_seg: Number of segments want to divide, can be 1 or larger than 1, but should be integer.
%     option: 'half'|'h' or 'whole'|'w', means you want a full frequency domain output or just half frequency domain output.
%
% Output:
%     freq: Output frequency,f = n*dt/N_segment, a N_segment*1 matrix, where n = [0:1:N_segment-1]';
%     apsd: Mean Auto-spectral Density(mean of M_segments).
%     
%   Author: Zelun Wu,
%   Ph.D. student of Physical Oceanography
%   University of Delaware, Xiamen University
%   zelunwu@udel.edu, zelunwu@stu.xmu.edu.cn
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%% Error test
switch option
    case {'half','h','whole','w'}
        ;
    otherwise
        error('Unknown "option" choice');
end

if rem(M_seg,1)>1e-10
    error('M_seg should be INTEGER');
elseif M_seg<1
    error('M_seg should be >= 1');
end
%-------------------------------------------------------------------------%
%% Calculation
N = length(ts); % Time series length.
N_seg = floor(N/M_seg); % Time series segment length. floor(x): rounddown function. E.g., floor(2.1) = 2
ts_seg_matrix = reshape(ts(1:N_seg*M_seg),[N_seg,M_seg]); % Reshape the first N_seg*M_seg(tails are abandoned) numbers into a N_seg*M_seg marix. Doing this step to eliminate loop computation.
T_seg = N_seg*dt; % Time length of segemet.
%-------------------------------------------------------------------------%
%% Discrete Fourier Transform
[Xf_seg_matrix,freq] = o_dft(ts_seg_matrix,dt,option); % Discrete Fourier Transform of ts segment matrix. if M_seg>1, the outour Xf will be a N_seg*M_seg matrix.
%-------------------------------------------------------------------------%
%% Auto Spectrum Density. APSD = 1/T*|DFT(ts)|^2
apsd_seg_matrix = 1./T_seg*abs(Xf_seg_matrix).^2; % abs(Xf) = real(Xf).^2 + imag(Xf).^2, where Xf is complex number.
apsd = mean(apsd_seg_matrix,2); % Mean APSD.
%-------------------------------------------------------------------------%
end