function [Xf,freq] = o_dft(ts,dt,option)
%%
%  Xf = O_DFT(ts,dt)
%  Discrete Fourier Transform of time series ts with sample rate dt (time domain)
%  Input:
%     ts: time series, should be a N(sample_points)*1 sequence, or N(sample_points) * M_segment matrix. M_segment should be >= 1.
%     dt: sample rate
%     option: 'h' or 'w', stands for half frequency domain or the whole frequency domain
%
%   Output:
%     freq: frequency
%     Xf: Fourier Transform of ts in frequency domain. The same size as ts, or N(sample_points)/2 * M_segment(if option = 'half');
%     
%   Author: Zelun Wu,
%   Ph.D. student of Physical Oceanography
%   University of Delaware, Xiamen University
%   zelunwu@udel.edu, zelunwu@stu.xmu.edu.cn

%-------------------------------------------------------------------------%
%% Error test
switch option
    case {'half','h','whole','w'}
        ;
    otherwise
        error('Unknown "option" choice');
end

[sz1,sz2] = size(ts); %sz1 should be Time series (can be segment) length, and sz2 should be M_segment.
if sz1<sz2
    warning('The number of columns(should be number of segments) is larger than number of rows(should be length of segemts)');
end
%-------------------------------------------------------------------------%
%% Parameters
N = sz1; M = sz2; % Number of ts points
T = N * dt; % Sample length in temporal dimention
n = [0:N-1]'; % n = [0,1,2,...,N-1]', is a N by 1 matrix.
% calculate frequency: f(n) = n/(N*dt) = n/T, with n = 0,1,2, ... N-1      (1)
freq = n/T; % frequency
%-------------------------------------------------------------------------%
%% Calculate Xf: 
% Xf(f = n/T) = dt * sum_{k=0}^{N-1){x(k*dt) exp(-i*2pi(n*k)/N)}       (2),
% in the form of matrix computation: 
%             (     [ 0 ]                        ) [x_seg1(t=0*dt)     x_seg2(t=0*dt)     ...  x_segM(t=0*dt)    ]
%             (     [ 1 ]                        ) [x_seg1(t=1*dt)     x_seg2(t=1*dt)     ...  x_segM(t=1*dt)    ]
%             (     [ 2 ]                        ) [x_seg1(t=2*dt)     x_seg2(t=2*dt)     ...  x_segM(t=2*dt)    ]
% Xf = dt* exp(2*pi*[...]*[0,1,2,3,...,n,...,N-1])*[     ...              ...             ...        ...         ]
%             (     [ k ]                        ) [x_seg1(t=k*dt)     x_seg2(t=k*dt)     ...  x_segM(t=1*dt)    ]
%             (     [...]                        ) [     ...              ...             ...        ...         ]
%             (     [N-1]                        ) [x_seg1(t=(N-1)*dt) x_seg2(t=(N-1)*dt) ...  x_segM(t=(N-1)*dt)]
                           
k = [0:N-1]; % k = [0,1,2,...,N-1], is a 1 by N matrix.
arg = -i*2*pi/N*(n*k); % arg = (2*pi/N)C * n * k, is a N by N matrix, where the nth row vector = [n*0, n*1, n*2, ..., n*k, ..., n*N-1].
Xf = dt*(exp(arg)*ts); % Eq.(2) in a matrix form, that ts is a N by 1 matrix,  exp(-i*arg) has the same size with arg is a N by N matrix.
%-------------------------------------------------------------------------%
switch option
    case {'half','h',0.5}
        in = ceil(N/2); % ceil(): round-up function. E.g., ceil(5.1) = 6
        freq = freq(1:in,:);
        Xf = Xf(1:in,:);
    case {'whole','w',1}
        ; % Do noting.
end
%-------------------------------------------------------------------------%
end