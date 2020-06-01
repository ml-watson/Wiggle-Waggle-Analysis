%% After Analysing some flipping data with the ANALYSIS_FULL_WAGGLE script,
%this code will calculate and graph the optical torque function.
%Lachlan J. Gibson, lachlan.gibson@uqconnect.edu.au
%Last modified 08/02/2013


%% Input Variables
flip=transpose(real(Mphi2)); %Dependent Variable (angle)
flipunc=transpose(SMphi1)/(0.5*NumFlips)^0.5; %Dependent variable uncertainty
ttt=transpose(time); %Independent Variable (time)
nlin=51; %The number of points included in each fitted segment (odd number)
minangle=0.04; %The smallest angle calculated.
cutspace=20; %spacing between points on the final graph


%% Calculation of the gradient as a function of angle
%by linearly fitting sections of the data.
%Truncate the data at the minimum angle
[~,flipend]=min(abs(flip-minangle)); %Find the index of the minimum angle
flip=flip(1:flipend); %truncate the dependent variable
ttt=ttt(1:flipend); %truncate the independent variable
flipunc=flipunc(1:flipend); %truncate the error in independent variable

%preallocate memory
phi=flip((nlin+1)/2:end-(nlin-1)/2); %The angle variable
phidot=zeros(1,length(phi)); %The gradient
phidotunc=zeros(1,length(phi)); %Error in gradient

%This for loop calculates the gradient and error of the angle over time
for ii=1:length(phi)
    %Take the fit of segment
    fitobj=fit(ttt(ii:ii+nlin-1),flip(ii:ii+nlin-1),'poly1');
    %Extract the values from the fit
    coeff=coeffvalues(fitobj); %gradient and y intercept
    conf=confint(fitobj); %95% confidence intervals
    phidot(ii)=coeff(1); %The gradient by itself
    phidotunc(ii)=0.5*abs(conf(1,1)-conf(2,1)); %Error: half conf interval
end

%Calculation of gradient using the 2 end points only
phidot2=1e4*(flip(nlin:end)-flip(1:end-nlin+1))/(nlin-1); %gradient
phidot2unc=1e4*(flipunc(nlin:end).^2+...
    flipunc(1:end-nlin+1).^2).^0.5/(nlin-1); %uncertainty


%% Calculation of the decay constant (chi/8pia^2eta)
%The angular time derivative as a function of angle is fit with 0.5ksin(2x)
sinfit=fit(phi,transpose(phidot),fittype(@(K,x) -0.5*K*sin(2*x)),...
    'StartPoint',50); %fit
%Extract the decay constant with error from the fit
conf=confint(sinfit); %95% confidence bounds
Kunc=0.5*abs(conf(1)-conf(2)); %error is half confidence interval
K=coeffvalues(sinfit); %Decay Constant


%% The optical torque function is plotted
%The torque function is proportional to the angular time derivative
%phidot=-K*T(phi)

%plot new theory (sinusoidal torque)
plot((0:0.001:0.92),0.5*sin(2*(0:0.001:0.92)),'k')
hold on
%plot old theory (linear torque)
plot([0 0.92],[0,0.92], 'k--')
%plot data calculated from the two end points
errorbar(Cutter(phi,cutspace,1),Cutter(phidot2,cutspace,1)/-K,...
    Cutter(phidot2unc,cutspace,1)/K,'bo','MarkerSize',4.0)
%plot data calculated from linear fits
% errorbar(Cutter(phi,cutspace,1),Cutter(phidot,cutspace,1)/-K,...
%     Cutter(phidotunc,cutspace,1)/K,'ro','MarkerSize',4.0)
%Only every 'cutspace' point is plotted for the sake of clarity
xlabel('Angle ($\phi$, rad)','Interpreter','LaTex','FontSize',16)
ylabel('Optical Torque Function T($\phi$)',...
    'Interpreter','LaTex','FontSize',16)
legend('0.5sin(2\phi)','\phi','Experimental Data','Location','NorthWest')
title('Nonlinear Optical Torque','Interpreter','LaTex','FontSize',16)








