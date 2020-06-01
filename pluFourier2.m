function [omega, F] = pluFourier2(t,y,alpha)
%A piecewise linear unilateral Fourier transform (algorithm based on
%one by Tassieri). Data begin at time zero and extend sufficiently long
%that discretisation error is not too large.
%Inputs must be column vectors.
%
%This variant of pluFourier.m is designed to replace the linear
%extrapolation to infinity with an exponential decay, thereby ensuring that
%the physical limits of the autocorrelation function are satisfied.
%The parameter alpha describes the exponential decay
%               A = exp(-alpha*t)
%used to extrapolate to infinity.

%% Create a frequency grid.
t_span = max(t)-min(t);
del_t = t(2)-t(1);
%Create a frequency grid from the lowest harmonic to the Nyquist frequency.
N = 1e3; %Number of output points.
nu_min = 1./(2*t_span); %Twice the period of the sampling time.
nu_max = 1./(2*pi*del_t); %Half the sampling frequency.
nu = logspace(log10(nu_min),log10(nu_max),N);
omega = 2*pi*nu; %Convert to angular frequencies.

%% Initialise loops.
F = zeros(1,N); %Blank storage.
t(1) = []; %Remove the time zero point to get the correct phases.
 %The gradients between each consecutive pair of data points.
grads = diff(y)./del_t;
%Begin loops.
for jj = 1:N; %For every data point...
    out = grads.*exp(-1i*omega(jj).*t);
    F(jj) = (exp(1i*omega(jj)*del_t)-1)*sum(out);
end

%% Add the boundary terms.
F = F./((1i.*omega).^2) + ((y(1)-y(length(y))*exp(-1i*omega*max(t)))./(1i.*omega)) +...
    y(length(y))*exp(-1i*omega*max(t))./(alpha+1i.*omega);
end