function [freq,cpsd_xy_mean,coherence_xy,phase_x_mean,phase_y_mean] = csa(x,y,dt,M_seg,option)
% [freq,cpsd_xy_mean,coherence_xy,phase_x_mean,phase_y_mean] = csa(x,y,dt,M_seg,option)
% Cross Spectrum Analysis
% Based on the Discrete Fourier Transform function: dft.m.
% %%
% freq(n+1) = n*dt/N_segment, n = 0,1,2,...,N_segment-1
% cpsd = 1/Time_segment*conj(Xf)*Yf, where Xf and Yf is the Fourier Transform of x and y.
% If M_seg>1, then cpsd is the mean apsd of M_segments divided from x and y.
% %%
% Input:
%     x : Input time series 1. shoule be an N*1 or 1*N sequence.
%     y : Input time series 2. shoule be an N*1 or 1*N sequence.
%     dt: Sampling rate in time domain. E.g., 1 hour.
%     M_seg: Number of segments want to divide, can be 1 or larger than 1, but should be integer.
%     option: 'half'|'h' or 'whole'|'w', means you want a full frequency domain output or just half frequency domain output.
% Output:
%     freq: Output frequency,f = n*dt/N_segment, a N_segment*1 matrix, where n = [0:1:N_segment-1]';
%     cpsd_mean: Mean Cross Spectral Density(mean of M_segments).
%     coherence_xy: Mean Manitude Square Coherence of x and y
%     phase_x_mean: Mean phase of time series x.
%     phase_y_mean: Mean phase of time series y.
%     
%   Author: Zelun Wu,
%   Ph.D. student of Physical Oceanography
%   University of Delaware, Xiamen University
%   zelunwu@udel.edu, zelunwu@stu.xmu.edu.cn

%-------------------------------------------------------------------------%
%% Error Checking
N1 = length(x); N2 = length(y);
if N1 ~= N2;error('The length of two time series should be the same');end;

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
%%
N = N1;
T = N*dt;

N_seg = floor(N/M_seg);
T_seg = N_seg*dt;

x_seg_matrix = reshape(x(1:N_seg*M_seg),[N_seg,M_seg]);
y_seg_matrix = reshape(y(1:N_seg*M_seg),[N_seg,M_seg]);
%-------------------------------------------------------------------------%
%% cpsd
[Xf,freq] = dft(x_seg_matrix,dt,option);
[Yf,freq] = dft(y_seg_matrix,dt,option);

cpsd_xy_seg = 1/T_seg * conj(Xf).*Yf;
cpsd_xy_mean = mean(cpsd_xy_seg,2);
% cpsd_xy_mean = nanmean(cpsd_xy_seg_matrix,2);
% cpsd_xy_mean = 1/T_seg * diag(Yf_seg_matrix * Xf_seg_matrix')/M_seg;
%-------------------------------------------------------------------------%
%% coherence
% apsd_xx_mean = 1/T_seg * diag(Xf_seg_matrix*Xf_seg_matrix')/M_seg;
% apsd_yy_mean = 1/T_seg * diag(Yf_seg_matrix*Yf_seg_matrix')/M_seg;
apsd_xx_seg = 1./T_seg*abs(Xf).^2;
apsd_yy_seg = 1./T_seg*abs(Yf).^2;

apsd_xx_mean = mean(apsd_xx_seg,2);
apsd_yy_mean = mean(apsd_yy_seg,2);
coherence_xy = abs(cpsd_xy_mean).^2./(apsd_xx_mean.*apsd_yy_mean);
%-------------------------------------------------------------------------%
%% phase
phase_x_seg_matrix = atan(imag(Xf)./real(Xf));
phase_y_seg_matrix = atan(imag(Yf)./real(Yf));
phase_x_mean = nanmean(phase_x_seg_matrix,2);
phase_y_mean = nanmean(phase_y_seg_matrix,2);
%-------------------------------------------------------------------------%
end