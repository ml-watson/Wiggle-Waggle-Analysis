%% TDMS extractor.
% Uses the .tdms to .mat converter (convertTDMS) to create a .mat file
% containing appropriately placed and referenced data. Run this once to
% import data into the workspace with the correct formatting for use with
% ANALYSIS_FULL. After the initial run it is possible to simply load
% directly from .mat without this program.
%
%Data must be collected as two .tdms files; the TDMS convertor will not
%handle appended files!
%The first file contains measurements of the rotating sphere, for
%calibration of the angle from the HeNe polarisation difference channel.
%The second should contain the data for analysis (which should NOT be
%collected immediately after the calibration data: some time is required
%for the sphere and fluid to come to steady state).

%% Inputs.
% The program requires the filename from the working directory
fpath = '110316\water_m30_s40_v1';    %The data set itself.
fpath_waggle = [fpath '_waggle']; %The separate waggle file.
fpath_brown1 = [fpath '_brown1']; %The separate brownian 1 file.
fpath_brown2 = [fpath '_brown2']; %The separate brownian 2 file.
fpath_calib = [fpath '_calib']; %The separate calibration file.
%Also provide the apparent (on-screen) radius [cm] and uncertainty in the radius
%along with an estimate of the trapping power and ambient temperature.
a = 1.5; %Radius (micrometers, cm on screen).
del_a = 0.4/2; %Absolute uncertainty.
P_trap = [30 40]; %Estimated trapping power [mW] for slave (S) and master (M).
                    %(2.1W output, AOM calibration by Emmanuel Brousse).
T_amb = 273.15+23; %Lab temperature [K].
del_T_amb = 3; %Absolute uncertainty.
flip_time=10000; %Time between flips [ms] (flip once then flip back counts as two flips)
%Add a description of the experiment.
desc = 'Shu''s measurement in water on 11/03/2016';

%% The waggle data must be extracted from the .tdms format.
%This will write a .mat file with the name given in 'fpath'.
%IT WILL NOT GIVE WARNING OF AN OVERWRITE!
path_dat = [fpath_waggle '.tdms'];
convertTDMS(1,char(path_dat));

%% We now manipulate the structure that was produced.
%Load the .mat file.
fpath_waggle = [fpath_waggle '.mat'];
load(char(fpath_waggle));

% Form all of these data into a new structure.
% LVE --- Parameters --- Radius, Temp, Power, NAPAF truncation index
%     --- WAGGLE     --- t, HeNe_diff, HeNe_sum, IR.
LVE.Param.Radius = 2*5.0025e-7*[a del_a]; %Radius in microns.
LVE.Param.Temperature = [T_amb del_T_amb];
LVE.Param.Power = P_trap;
LVE.Param.NAPAFTrunc = []; %To be set manually later if desired.
LVE.Param.Description = desc;
LVE.Param.Fliptime= flip_time;
%Use an implicit timeseries, since the converter does not handle time
%correctly when implicitly defined in the .tdms.
LVE.WAGGLE.t = ConvertedData.Data.MeasuredData(5).Property(3).Value.*...
    (0:(length(ConvertedData.Data.MeasuredData(5).Data)-1))';
LVE.WAGGLE.HeNe_sum = ConvertedData.Data.MeasuredData(5).Data+ConvertedData.Data.MeasuredData(6).Data;
LVE.WAGGLE.HeNe_diff = ConvertedData.Data.MeasuredData(5).Data-ConvertedData.Data.MeasuredData(6).Data;
%LVE.WAGGLE.IR = ConvertedData.Data.MeasuredData(7).Data+ConvertedData.Data.MeasuredData(8).Data;
LVE.WAGGLE.IR = ConvertedData.Data.MeasuredData(4).Data;

clear ConvertedData

fprintf('\n Processing brownian data.')
    fprintf('\n')

%% The brownian data must be extracted from the .tdms format.
%This will write a .mat file with the name given in 'fpath'.
%IT WILL NOT GIVE WARNING OF AN OVERWRITE!
path_brown1 = [fpath_brown1 '.tdms'];
convertTDMS(1,char(path_brown1));
    
%Load the .mat file.
fpath_brown1 = [fpath_brown1 '.mat'];
load(char(fpath_brown1));

% Store all of these data into the previously made structure.
%LVE --- BROWN --- Data  --- t, HeNe_diff
%Use an implicit timeseries, since the converter does not handle time
%correctly when implicitly defined in the .tdms.
LVE.BROWN1.t = ConvertedData.Data.MeasuredData(5).Property(3).Value.*...
    (0:(length(ConvertedData.Data.MeasuredData(5).Data)-1))';
LVE.BROWN1.HeNe_sum = ConvertedData.Data.MeasuredData(5).Data+ConvertedData.Data.MeasuredData(6).Data;
LVE.BROWN1.HeNe_diff = ConvertedData.Data.MeasuredData(5).Data-ConvertedData.Data.MeasuredData(6).Data;
%LVE.BROWN1.HeNe_diff = LVE.BROWN1.HeNe_diff-mean(LVE.BROWN1.HeNe_diff);
%LVE.BROWN1.IR = ConvertedData.Data.MeasuredData(7).Data+ConvertedData.Data.MeasuredData(8).Data;
LVE.BROWN1.IR = ConvertedData.Data.MeasuredData(4).Data;

clear ConvertedData

fprintf('\n Processing calibration data.')
    fprintf('\n')

%% The brownian data must be extracted from the .tdms format.
%This will write a .mat file with the name given in 'fpath'.
%IT WILL NOT GIVE WARNING OF AN OVERWRITE!
path_brown2 = [fpath_brown2 '.tdms'];
convertTDMS(1,char(path_brown2));
    
%Load the .mat file.
fpath_brown2 = [fpath_brown2 '.mat'];
load(char(fpath_brown2));

% Store all of these data into the previously made structure.
%LVE --- BROWN --- Data  --- t, HeNe_diff
%Use an implicit timeseries, since the converter does not handle time
%correctly when implicitly defined in the .tdms.
LVE.BROWN2.t = ConvertedData.Data.MeasuredData(5).Property(3).Value.*...
    (0:(length(ConvertedData.Data.MeasuredData(5).Data)-1))';
LVE.BROWN2.HeNe_sum = ConvertedData.Data.MeasuredData(5).Data+ConvertedData.Data.MeasuredData(6).Data;
LVE.BROWN2.HeNe_diff = ConvertedData.Data.MeasuredData(5).Data-ConvertedData.Data.MeasuredData(6).Data;
%LVE.BROWN2.HeNe_diff = LVE.BROWN2.HeNe_diff-mean(LVE.BROWN2.HeNe_diff);
%LVE.BROWN2.IR = ConvertedData.Data.MeasuredData(7).Data+ConvertedData.Data.MeasuredData(8).Data;
LVE.BROWN2.IR = ConvertedData.Data.MeasuredData(4).Data;

clear ConvertedData

fprintf('\n Processing calibration data.')
    fprintf('\n')

%% The calibration data must be extracted from the .tdms format.
%This will write a .mat file with the name given in 'fpath'.
%IT WILL NOT GIVE WARNING OF AN OVERWRITE!
path_calib = [fpath_calib '.tdms'];
convertTDMS(1,char(path_calib));
    
%Load the .mat file.
fpath_calib = [fpath_calib '.mat'];
load(char(fpath_calib));

% Store necessary calibration data into the previously made structure.
%LVE --- CALIB --- Calib_diff
LVE.CALIB.Calib_diff = ConvertedData.Data.MeasuredData(5).Data-ConvertedData.Data.MeasuredData(6).Data;

%% Write the final structure to disc, clearing the previous two.
save(fpath,'LVE')
delete(fpath_calib)
delete(fpath_brown1)
delete(fpath_brown2)
delete(fpath_waggle)
clear all
fprintf('\n Restructuring complete (overwritten previous .mat).')
    fprintf('\n')