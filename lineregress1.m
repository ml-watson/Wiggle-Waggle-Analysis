function param = lineregress1(x,y)
%Linear regression of equally weighted data.
%Takes input x and y.

D = sum((x - mean(x)).^2); %Find D.
m = (1/D)*sum((x - mean(x)).*y); %Find m.
c = mean(y) - m*mean(x); %Find c.

%Calculate the squared residuals.
d_sq = (y - m*x-c).^2;

%Find uncertainties.
delta_m2 = (1/(D*(length(x)-2)))*sum(d_sq);
    delta_m = sqrt(delta_m2); %Uncertainty in slope.
delta_c = sqrt((1/length(x) + (mean(x)^2)/D)*((1/(length(x)-2)))*sum(d_sq)); %Uncertainty in intercept.

param = [m delta_m; c delta_c]; %Output = [ slope, uncertainty in slope; intercept, uncertainty in intercept].
end