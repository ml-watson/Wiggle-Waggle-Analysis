%% Waggle Analysis Program
% Determines G* by analysing the rotational dynamics of a vaterite flipping
% between two angularly offset traps.
% Originally Written by Lachlan J. Gibson, lachlan.gibson@uqconnect.edu.au
% Updates written by Mark Watson, mark.watson@uq.edu.au
% Last modified 1st June 2020 by Mark Watson


%Input Data
if exist('LVE','var') == 0
    disp('Loading Data')
    load('R:\SMP\Research\OMG\Lachlan Gibson\Matlab\Wiggle_Waggle\Waggle Scripts and Functions\110316\water_m30_s40_v1.mat')
end


%% Options
BM=0; % Set to 1 to analyse the Brownian Motion using autocorrelation
calibfail=0; % Set to 1 to correct calibration fail
varprogress=0; % Set to 1 to see how the variance changes throughout time
no_normalise=1; % Set to 1 to see averaging of flips without normalisation
plotting=1; % Set to 1 to plot complex shear modulus


%% Angle calibration
% The long axis of the transmitted elliptically polarised HeNe light is 
% parallel to the optic axis of the vaterite probe particle.
% The orientation of the polarisation ellipse determines the relative
% powers of perpendicular linear components. The difference in these powers
% is related sinusoidally to the azimuthal angle of the vaterite's optic
% axis: dV = A sin(2 phi) + B. The coefficients A and B are found from the
% maximum and minimum voltage difference.

dV=LVE.CALIB.Calib_diff; % The voltage difference during calibration
A=(max(dV)-min(dV))/2; % The amplitude is half the range
B=(max(dV)+min(dV))/2; % The offset is the average value (middle of range)

% The probe orientations can now be calculated using the calibration
phibrown1=(0.5*asin((LVE.BROWN1.HeNe_diff-B)/A));
phibrown2=(0.5*asin((LVE.BROWN2.HeNe_diff-B)/A));
phiwaggle=(0.5*asin((LVE.WAGGLE.HeNe_diff-B)/A));


%clear unnecessary variables
clear A B dV


%% Parameters and Constants
% The radius [m].
a=LVE.Param.Radius;
% Find the trap stiffness by equipartition.
k=1.38065*10^-23; %Boltzmann constant [J/K].
chi01=k*LVE.Param.Temperature(1)/var(phibrown1); %Stiffness of AOM 1
chi02=k*LVE.Param.Temperature(1)/var(phibrown2); %Stiffness of AOM 2


%% Find the precise start time of each flip using the IR data
%Prepare the IR data for a histogram by removing most of the excess data
IR=LVE.WAGGLE.IR;
IR=IR(IR>(2*mean(IR)-max(IR))/(4/3));

% Set the mean to zero so the AOMs data is distinguishable by sign
IR=IR-mean(IR);
%IR=IR(1:2430000)
% Duplicate with 1 timestep lag.
IR2 = IR(1:(length(IR)-1));
IR1 = IR(2: length(IR));

% Indicies when the AOMs switch is given when IR*IR2 is negative
indi=find(sign(IR1.*IR2)==-1)+1;

% clear unnecessary varibles
clear IR2 IR IR1


%% Analysis of Brownian Motion (Wiggle Method) for Comparison
% If the user sets BM=1 then the Brownian Motion from each trap is analysed
if BM==1
    % The time vectors for brownian motion
    t1=LVE.BROWN1.t; 
    t2=LVE.BROWN2.t;
    
    % The normalised autocorrelation functions
    A1=autocorr(phibrown1,length(t1)-1);
    A2=autocorr(phibrown2,length(t2)-1);
    
    % NAPAF 1 is plotted so the user can pick a truncation time
    figure(3),
    screen_size = get(0, 'ScreenSize');
    f3=figure(3);
    set(f3, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    semilogx(A1),
    xlabel('Time (1e-4 s)','FontSize',24),
    ylabel('NAPAF 1','FontSize',24),
    title('Please provide the truncation point')
    [trunc1,~] = ginput(1); %User indicates the truncation point.
    close 3
    
    % NAPAF 2 is plotted so the user can pick a truncation time
    figure(4),
    screen_size = get(0, 'ScreenSize');
    f4=figure(4);
    set(f4, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    semilogx(A2),
    xlabel('Time (1e-4 s)','FontSize',24),
    ylabel('NAPAF 2','FontSize',24),
    title('Please provide the truncation point')
    [trunc2,~] = ginput(1); %User indicates the truncation point.
    close 4
    
    
    % Round the tuncation indices to integers
    trunc1=floor(trunc1); trunc2=floor(trunc2);
    
    % Truncate the NAPAFs and times
    A1=A1(1:trunc1); A2=A2(1:trunc2);
    t1=t1(1:trunc1); t2=t2(1:trunc2);
    
    % Linear fit the truncation foint to estimate the exponential decay
    param=lineregress1(t1((trunc1-100):trunc1),log(A1((trunc1-100):trunc1)));
    alpha_est=abs(param(1,1));
    [omega1,B1]=pluFourier2(t1,A1,alpha_est);
    
    % Also linear fit the second data set
    param=lineregress1(t2((trunc2-100):trunc2),log(A2((trunc2-100):trunc2)));
    alpha_est=abs(param(1,1));
    [omega2,B2]=pluFourier2(t2,A2,alpha_est);
    
    % Calculate the complex shear modulus from the Fourier transforms
    G_star01=-(1/(8*pi*(a(1)^3)))*(chi01*omega1.*B1./(1i+omega1.*B1));
    G_star02=-(1/(8*pi*(a(1)^3)))*(chi02*omega2.*B2./(1i+omega2.*B2));
    
end



%% Segmate each flip based on the flip times
%For the sake of averaging, each flip is truncated at the time of the
%shortest flip
NumFlips=length(indi)-1; %The total number of flips
fliplength=min(diff(indi)); %The length of the shortest flip

%Preallocate memory for variables used in the for loop
flipphi=zeros(NumFlips,fliplength); %A 2D matrix ready to contain each flip
point=zeros(1,NumFlips);

% The direction of each flip is established
phiwaggle=phiwaggle-mean(phiwaggle);
%phiwaggle=phiwaggle(1:2430000);
flipdirection=sign(phiwaggle(indi)); % - means up and + means down

for ii=1:NumFlips % the loop cycles for each flip
    
    % import the data into the 2D marix flipphi
    flipphi(ii,:)=phiwaggle(indi(ii):indi(ii)+fliplength-1);
    
    % set all flips to the same direction (every second one is backwards)
    flipphi(ii,:)=flipphi(ii,:)*flipdirection(ii);
    
    % An initial estimate of when the flip finishes is given by the first
    % point when the flip reaches the mode of its angular distribution.
    [H,Hphi]=ksdensity(flipphi(ii,:)); %the kernel density distribution
    point(ii)=find(flipphi(ii,:)<Hphi(H==max(H)),1); %find the point at the mode
end

% clear unnecessary varibles
clear indi H Hphi flipdirection


%% CALIBRATION-FAIL-CORRECTION
%If the vaterite flips past the calibration maximum this will correct it
if calibfail==1
    turningpoint=zeros(1,NumFlips);
    for ii=1:NumFlips %the loop cycles for each flip
        turningpoint(ii)=find(flipphi(ii,:)==max(flipphi(ii,:)),1);
        flipphi(ii,1:turningpoint(ii))=-flipphi(ii,1:turningpoint(ii))+2*max(flipphi(ii,:));
    end
end


%% Analysis of variance using the motion between flips
% Lower and upper bounds of the trap stiffness can be found by stitching
% together the BM at the end of each flip. This will over estimate the
% variance. Mean aligning each section, however, gives a lower estimate of
% the variance.
%Preallocate memory for variables used in the for loops
phibrown11=zeros(1,sum(fliplength-point(1:2:end)+1));
phibrown12=zeros(1,sum(fliplength-point(2:2:end)+1));
phibrown13=zeros(1,sum(fliplength-point(1:2:end)+1));
phibrown14=zeros(1,sum(fliplength-point(2:2:end)+1));
n=1e2; % lag time interval
var_phi11=zeros(NumFlips,round(fliplength/n));
var_phi=zeros(1,NumFlips);

%The BM data from the 'flipphi' matrix is imported into another vector(s)
indc=1;
for ii=1:2:NumFlips % odd flips
    indl=fliplength-point(ii)+1; % length of single flip BM data
    % Import the mean shifted flip data
    phibrown11(indc:indc+indl-1)=flipphi(ii,point(ii):fliplength)...
        -mean(flipphi(ii,point(ii):fliplength));
    % Import the flip data (without mean shifting)
    phibrown13(indc:indc+indl-1)=flipphi(ii,point(ii):fliplength);
    indc=indc+indl; % indc keeps track of the index
end
indc=1;
for ii=2:2:NumFlips % even flips
    indl=fliplength-point(ii)+1; % length of BM flip data
    % Import the mean shifted flip data
    phibrown12(indc:indc+indl-1)=flipphi(ii,point(ii):fliplength)...
        -mean(flipphi(ii,point(ii):fliplength));
    % Import the flip data (without mean shifting)
    phibrown14(indc:indc+indl-1)=flipphi(ii,point(ii):fliplength);
    indc=indc+indl; % indc keeps track of the index
end

if varprogress==1
    % Extra variance analysis (not necessary)
    for ii=1:NumFlips % the loop cycles for each flip
        % variance at differnt lag times
        for jj=n:n:fliplength-point(ii)
            var_phi11(ii,jj/n)=var(flipphi(ii,point(ii):point(ii)+jj));
        end
        
        %find the individual variances for each flip (after 'point(ii)')
        var_phi(ii)=var(flipphi(ii,point(ii):fliplength));
    end
end


%The trap stiffness from the flipping data
chi11=k*LVE.Param.Temperature(1)/var(phibrown11);
chi12=k*LVE.Param.Temperature(1)/var(phibrown12);
chi13=k*LVE.Param.Temperature(1)/var(phibrown13);
chi14=k*LVE.Param.Temperature(1)/var(phibrown14);

%The actual trap stiffness is about half way between the two bounds
chi1=(chi11+chi13)/2;
chi2=(chi12+chi14)/2;

%error in trap stiffness is about half the range of the bounds
chi1unc=(chi11-chi13)/2;
chi2unc=(chi12-chi14)/2;

% clear unnecessary varibles
clear chi11 chi12 chi13 chi14 indc indl


%% More flip Preparation
flipheight=zeros(1,NumFlips);
for ii=1:NumFlips %the loop cycles for each flip
    %Set the equilibrium (the mean after the 'point') to zero
    flipphi(ii,:)=flipphi(ii,:)-mean(flipphi(ii,point(ii):fliplength));
    
    %The maximum angle of each flip
    flipheight(ii)=max(flipphi(ii,:));
end


%% Averaging the Flips Without Normalisation (only valid for purely viscous fluids)
if no_normalise==1
    fstart1=min(flipheight(1:2:end)); %The max common height acheived by odd flips
    fstart2=min(flipheight(2:2:end)); %The max common height acheived by even flips
    flipstart=zeros(1,NumFlips); %Variable which stores the index of common height
    flipphicut=zeros(NumFlips,fliplength); %New flip variable (similar to flipphi)
    
    %Truncate the odd flips at the common height
    for ii=1:2:NumFlips
        %find the index of the point closest to the hight
        [~,flipstart(ii)]=min(abs(fstart1-flipphi(ii,1:fliplength)));
        %store the flip data after truncating in flipphicut
        flipphicut(ii,1:fliplength-flipstart(ii)+1)=flipphi(ii,flipstart(ii):fliplength);
    end
    
    %Truncate the even flips at the common height
    for ii=2:2:NumFlips
        %find the index of the point closest to the hight
        [~,flipstart(ii)]=min(abs(fstart2-flipphi(ii,1:fliplength)));
        %store the flip data after truncating in flipphicut
        flipphicut(ii,1:fliplength-flipstart(ii)+1)=flipphi(ii,flipstart(ii):fliplength);
    end
    
    
    % Average the flips at the same starting angle
    NMphicut1=mean(flipphicut(1:2:end,:)); % odd mean
    STphicut1=std(flipphicut(1:2:end,:));  % odd stdv
    NMphicut2=mean(flipphicut(2:2:end,:)); % even mean
    STphicut2=std(flipphicut(2:2:end,:));  % even stdv
    
    %Average the flips at time zero for maximum angle range
    Mphi1=mean(flipphi(1:2:end,:)); % odd mean
    SMphi1=std(flipphi(1:2:end,:)); % odd stdv
    Mphi2=mean(flipphi(2:2:end,:)); % even mean
    SMphi2=std(flipphi(2:2:end,:)); % even stdv
    
    % clear unnecessary varibles
    clear flipheight fstart1 fstart2 flipstart flipphicut
end


%% Transformation, Normalisation and Averaging of Flips
for ii=1:NumFlips
    %trasformation to compensate for non quadratic trap potential
    flipphi(ii,1:fliplength)=tan(flipphi(ii,1:fliplength));
    
    %finding an approximate initial value (better than the first point)
    %using a linear regression of the first 20 points (code fromlineregress1)
    D = sum(((0:19) - mean((0:19))).^2); %Find D.
    m = (1/D)*sum(((0:19) - mean((0:19))).*flipphi(ii,1:20)); %Find m.
    c = mean(flipphi(ii,1:20)) - m*mean((0:19)); %Find c.
    %divide through by inital value
    flipphi(ii,1:fliplength)=flipphi(ii,1:fliplength)/c;
end


%find the mean and standard deviation of the transformed flips
NMTAP1=mean(flipphi(1:2:end,:)); % odd mean
STTAP1=std(flipphi(1:2:end,:));  % odd stdv
NMTAP2=mean(flipphi(2:2:end,:)); % even mean
STTAP2=std(flipphi(2:2:end,:));  % even stdv

% clear unnecessary varibles
clear D m c


%% Calculate G* from the transformed Normalised Mean Transformed Angular Position (NMTAP)
%Use the piecwise linear unilateral transform of Tassieri et al.
time=1e-4*(0:fliplength-1); %Time, assuming a time-step of 0.1ms
[omega, B1] = pluFourier(time,NMTAP1,[1 0]); %The transform of flip 1
[~, B2] = pluFourier(time,NMTAP2,[1 0]); %The transform of flip 2

%G* assuming low Reynolds number regime. For viscoelastic fluids the memory
%decay time should be much shorter than the flip decay time. (Shown to be
%valid for our experiment through simulations, see page 35 of LG's book 3)
G_star1 = -(1/(8*pi*(a(1)^3)))*(chi01*omega.*B1./(1i + omega.*B1));
G_star2 = -(1/(8*pi*(a(1)^3)))*(chi02*omega.*B2./(1i + omega.*B2));


% clear unnecessary varibles
%clear B1 B2 fliplength ii jj n


%% Plotting - viscoelasticity modulus
if plotting==1
    figure()
    loglog(omega,imag(G_star1),'k')
    hold on
    loglog(omega,real(G_star1),'r')
    loglog(omega,imag(G_star2),'k-.')
    loglog(omega,real(G_star2),'r-.')
    
    if BM==1
        loglog(omega1,imag(G_star01),'b')
        loglog(omega1,real(G_star01),'g')
        loglog(omega2,imag(G_star02),'b-.')
        loglog(omega2,real(G_star02),'g-.')
        legend('waggle_odd loss', 'waggle_odd storage','waggle_even loss','waggle_even storage','brown1 loss', 'brown1 storage','brown2 loss','brown2 storage','Location', 'Northwest') 
    else
        legend('waggle_odd loss', 'waggle_odd storage','waggle_even loss','waggle_even storage','Location', 'Northwest') 
    end
end

%% extra ploting
figure()
plot(phiwaggle)



