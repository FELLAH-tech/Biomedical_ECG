function [X, f, t] = stft(x,w,d,N_fft,Fs)
% This function computes the stft for m = [0, d, 2d, 3d...]
% This function outputs are:
% -> X, which is a matrix of n_fft lines and M columns
% M is the number of elements of m
% X(i,j) is the value of the spectrogram for time t(i) and frequency f(j)
% -> f, is a column vector of the frequencies (in Hz)
% -> t, is a row vector containing the times of the beginning of the windows
L = length(x);
N = length(w);
M = floor((L-N)/d)+1;
Y = zeros(N,M);
X = zeros(N_fft,M);
w = w(:); % make sure it's a column
c=1;
for i = 1 : M
    if i==255
        c=c+1;
    end
    if (i == 1)
        Y(:,1) = x(1:N);
    else
        if (N+(i-1)*d < L)
            Y(:,i) = x((i-1)*d:N+(i-1)*d-1);% input data frame
            z = w.*Y(:,i); % Apply window
            X(:,i) = fft(z,N_fft);

        else
            Y(:,M) = [x((M-1)*d:L-1) zeros(1,N-L+(M-1)*d-1)];
            z_m = w.*Y(:,M);
            X(:,M) = fft(z_m,N_fft);
        end
    end
end
t = (0:M-1)*d/Fs;
f1 = 0:1/N_fft:(N_fft-1)/N_fft;
f=f1*Fs;



