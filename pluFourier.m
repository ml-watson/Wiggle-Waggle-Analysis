function [omega, F] = pluFourier(t,y,lims)
%A piecewise linear unilateral Fourier transform (algorithm by Tassieri).
%Data begin at time zero and extend sufficiently long that discretisation
%error is not too large. The integrand must approach g0 at zero and a
%gradient of g_inf at infinite times. Inputs must be column vectors.

%% %Define inputs.
%Remove the first data point. This is to ensure that the minimum time is
%non-zero.
t(1) = [];
y(1) = [];

%% Options, set for the NAPAF.
g0 = lims(1); %Value of the integrand at time zero.
g_inf = lims(2); %Limit of the gradient as t-> infty.

maxtt = max(t);
mintt = min(t);
%Create the upper and lower time vectors for rapid evaluation of the loop.
t_upper = t;
t_lower = t;
t_upper(1) = [];
t_lower(length(t)) = [];

%% Initialise loops.
N = 1e3; %Number of output points.
G0 = (y(1)-g0)/t(1);
%Create a frequency grid from the lowest harmonic to the Nyquist frequency.

nu_min = 1./(2*(maxtt-mintt)); %Twice the period of the sampling time.
nu_max = 1./(2*pi*(t(3)-t(2))); %Half the sampling frequency.
nu = logspace(log10(nu_min),log10(nu_max),N);
%nu = linspace(nu_min, nu_max,N);
omega = 2*pi*nu;
%omega = nu; %FAKE!

F = zeros(1,N); %Blank storage.
%% %Begin loops.
for jj = 1:N; %For every data point...
    out = (diff(y)./diff(t)).*(exp(-t_lower*1i*omega(jj))-exp(-t_upper*1i*omega(jj)));
    F(jj) = -(omega(jj)^(-2))*(1i*omega(jj)*g0 + ...
            g_inf*exp(-1i*omega(jj)*maxtt) + ...
            sum(out) + G0*(1-exp(-1i*omega(jj)*t(1))));
end

end