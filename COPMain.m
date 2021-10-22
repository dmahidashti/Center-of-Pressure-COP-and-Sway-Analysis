
%
%%

clear all; 
close all; 

load('/Volumes/GoogleDrive/My Drive/University/BME705 Rehabilitation/Lab 3/Data/Part 1/ECEO_data3.mat');
load ('/Volumes/GoogleDrive/My Drive/University/BME705 Rehabilitation/Lab 3/Data/Part 2/sway_data3.mat');

fs=2000;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUT = 16 channels of raw forceplate data;
% OUTPUT: xCOPL, yCOPL, xCOPR, yCOPR, Fz_R, Fz_L
% Note: Rv_R and Rv_L are vertical forces of right and left forceplate 
% Repeat for esEO_fp
[xCOPL_eo, yCOPL_eo, xCOPR_eo, yCOPR_eo, Rv_R_eo, Fz_L_eo] = analyzeFP(qsEO_fp);
[xCOPL_ec, yCOPL_ec, xCOPR_ec, yCOPR_ec, Rv_R_ec, Fz_L_ec] = analyzeFP(qsEC_fp);

%%% Calcualte the COPnet 
% 1st: Define the global center of pressure coodinate & use these to compute COPnet
% Shift X-coordinates to the middle of the split force plate by +/- 12.5825cm
xCOPLnew_eo = xCOPL_eo - 12.5825;
xCOPRnew_eo = xCOPR_eo + 12.5825;

xCOPLnew_ec = xCOPL_ec - 12.5825;
xCOPRnew_ec = xCOPR_ec + 12.5825;

% Calculate COP for overall system using Winter et al., 2003 [1]: 
% Equation A from the lab manual - use the new values adjusted by the size
% of the pressure plate
xCOP_eo = (xCOPLnew_eo.*(Fz_L_eo./(Fz_L_eo+Rv_R_eo))) + (xCOPRnew_eo.*(Rv_R_eo./(Fz_L_eo+Rv_R_eo)));
yCOP_eo = (yCOPL_eo.*(Fz_L_eo./(Fz_L_eo+Rv_R_eo))) + (yCOPR_eo.*(Rv_R_eo./(Fz_L_eo+Rv_R_eo)));

xCOP_ec = (xCOPLnew_ec.*(Fz_L_ec./(Fz_L_ec+Rv_R_ec))) + (xCOPRnew_ec.*(Rv_R_ec./(Fz_L_ec+Rv_R_ec)));
yCOP_ec = (yCOPL_ec.*(Fz_L_ec./(Fz_L_ec+Rv_R_ec))) + (yCOPR_ec.*(Rv_R_ec./(Fz_L_ec+Rv_R_ec)));

%%% Filter ALL data 
% Use a Low-pass filter: a Butterworth is recommended
% define corner frequency and normalized cut-off frequency
fc=10;
Wn = fc/(fs/2);
[b, a] = butter(4, Wn,'low');

% Eyes Open
    % X Component
    xCOPLnew_eo_filter = filter(b, a, xCOPLnew_eo); %filtering xCOP left while eyes open
    xCOPRnew_eo_filter = filter(b, a, xCOPRnew_eo); %filtering xCOP right while eyes open
    xCOP_eo_filter = filter(b, a, xCOP_eo); %filtering xCOP net while eyes open
    % Y Component
    yCOPL_eo_filter = filter(b, a, yCOPL_eo); %filtering yCOP left while eyes open
    yCOPR_eo_filter = filter(b, a, yCOPR_eo); %filtering yCOP right while eyes open
    yCOP_eo_filter = filter(b, a, yCOP_eo); %filtering yCOP net while eyes open

% Eyes Closed
    % X Component
    xCOPLnew_ec_filter = filter(b, a, xCOPLnew_ec); %filtering xCOP left while eyes close
    xCOPRnew_ec_filter = filter(b, a, xCOPRnew_ec); %filtering xCOP right while eyes close
    xCOP_ec_filter = filter(b, a, xCOP_ec); %filtering xCOP net while eyes close
    % Y Component
    yCOPL_ec_filter = filter(b, a, yCOPL_ec); %filtering yCOP left while eyes close
    yCOPR_ec_filter = filter(b, a, yCOPR_ec); %filtering yCOP right while eyes close
    yCOP_ec_filter = filter(b, a, yCOP_ec); %filtering yCOP net while eyes close

%%% Plot net COP (the stabilograms)

figure;
plot(xCOP_eo_filter, yCOP_eo_filter);
xlabel('M/L');
ylabel('A/P');
title('Center of Presseure: Stabiliogram for eyes open');

figure;
plot(xCOP_ec_filter, yCOP_ec_filter);
xlabel('M/L');
ylabel('A/P');
title('Center of Presseure: Stabiliogram for eyes closed');

%%% Using the method described in Prieto el al., 1994 paper [2]

% Select 20 seconds of the data in the time frame from 40-60 seconds
T=20;
%Convert Time to Index
t1 = 40*fs;
t2 = 60*fs;
% Selecting Array Portion
Ap20c = yCOP_ec(t1:t2-1, 1);
Ml20c = xCOP_ec(t1:t2-1, 1);
Ap20o = yCOP_eo(t1:t2-1, 1);
Ml20o = xCOP_eo(t1:t2-1, 1);
% Downsample to 100 Hz, using built in function
fs_new = 100;
% Re-writing selected portions
Ap20c = downsample(Ap20c, T);
Ml20c = downsample(Ml20c, T);
Ap20o = downsample(Ap20o, T);
Ml20o = downsample(Ml20o, T);
% Zero-mean AP and ML : Equations 1 & 2(Prieto el al., 1994)
APc = Ap20c - mean(Ap20c);
MLc = Ml20c - mean(Ml20c);
RDc = sqrt(APc.^2 + MLc.^2);% Resultant Distance for closed eyes
APo = Ap20o - mean(Ap20o);
MLo = Ml20o - mean(Ml20o);
RDo = sqrt(APo.^2 + MLo.^2);% Resultant Distance for open eyes

% Total Excursion: Equations 8 and 9 (Prieto el al., 1994) (For loop required here)
TOTEXap = 0;
TOTEXml = 0;
TOTEX = 0;
% Total Excursion for eyes closed
TOTEXc = 0;
for i = 1: length(APc)-1
    TOTEXc = TOTEXc + sqrt((APc(i+1)-APc(i))^2 + (MLc(i+1)-MLc(i))^2);
end
TOTEXapc = 0;
for i = 1: length(APc)-1
    TOTEXapc = TOTEXapc + APc(i+1)-APc(i) ;
end
TOTEXmlc = 0;
for i = 1: length(MLc)-1
    TOTEXmlc = TOTEXmlc + MLc(i+1)-MLc(i) ;
end

% Total Excursion for eyes open
TOTEXo = 0;
for i = 1: length(APo)-1
    TOTEXo = TOTEXo + sqrt((APo(i+1)-APo(i))^2 + (MLo(i+1)-MLo(i))^2);
end
TOTEXapo = 0;
for i = 1: length(APo)-1
    TOTEXapo = TOTEXapo + APo(i+1)-APo(i) ;
end
TOTEXmlo = 0;
for i = 1: length(MLo)-1
    TOTEXmlo = TOTEXmlo + MLo(i+1)-MLo(i) ;
end

% Mean velocity: Equations 10 and 11 (Prieto el al., 1994)
MVELOc = TOTEXc*fs;
MVELOapc = TOTEXapc*fs;
MVELOmlc = TOTEXmlc*fs;

MVELOo = TOTEXo*fs;
MVELOapo = TOTEXapo*fs;
MVELOmlo = TOTEXmlo*fs;

%% TOTX Display
TOT_disp =['The Total Excursions are as follows:'];
disp(TOT_disp);
TOTAP_disp =['In ECEO Data 3, the AP Total Excursion for Eyes Open is ',num2str(TOTEXapo)];
disp(TOTAP_disp);
TOTML_disp =['In ECEO Data 3, the ML Total Excursion for Eyes Open is ',num2str(TOTEXmlo)];
disp(TOTML_disp);
TOTN_disp =['In ECEO Data 3, the Net Total Excursion for Eyes Open is ',num2str(TOTEXo)];
disp(TOTN_disp);
TOTAP_disp =['In ECEO Data 3, the AP Total Excursion for Eyes Closed is ',num2str(TOTEXapc)];
disp(TOTAP_disp);
TOTML_disp =['In ECEO Data 3, the ML Total Excursion for Eyes Closed is ',num2str(TOTEXmlc)];
disp(TOTML_disp);
TOTN_disp =['In ECEO Data 3, the Net Total Excursion for Eyes Closed is ',num2str(TOTEXc)];
disp(TOTN_disp);
%% Velocity Display
TOT_disp =['The Total Mean Velocities are as follows:'];
disp(TOT_disp);
TOTAP_disp =['In ECEO Data 3, the AP Mean Velocities for Eyes Open is ',num2str(MVELOapo)];
disp(TOTAP_disp);
TOTML_disp =['In ECEO Data 3, the ML Mean Velocities for Eyes Open is ',num2str(MVELOmlo)];
disp(TOTML_disp);
TOTN_disp =['In ECEO Data 3, the Net Mean Velocities for Eyes Open is ',num2str(MVELOo)];
disp(TOTN_disp);
TOTAP_disp =['In ECEO Data 3, the AP Mean Velocities for Eyes Closed is ',num2str(MVELOapc)];
disp(TOTAP_disp);
TOTML_disp =['In ECEO Data 3, the ML Mean Velocities for Eyes Closed is ',num2str(MVELOmlc)];
disp(TOTML_disp);
TOTN_disp =['In ECEO Data 3, the Net Mean Velocities for Eyes Closed is ',num2str(MVELOc)];
disp(TOTN_disp);

%% Part 2

%Redifining the variables
sol_f=sol_fast;
sol_s=sol_slow;
ta_f=ta_fast;
ta_s=ta_slow;
%% Force plate Analysis
% This seciton should follow the first part of the analysis

% call the analyzeFP function again (do not edit)a

[xCOPL, yCOPL, xCOPR, yCOPR, Rv_R, Fz_L] = analyzeFP(fp_fast);
[xCOPL2, yCOPL2, xCOPR2, yCOPR2, Rv_R2, Fz_L2] = analyzeFP(fp_slow);

%%% Calcualte the COPnet
% 1st: Define the global center of pressure coodinate & use these to compute COPnet
% Shift X-coordinates to the middle of the split force plate by +/- 12.5825cm
xCOPLnew = xCOPL - 12.5825;
xCOPRnew = xCOPR + 12.5825;

xCOPLnew2 = xCOPL2 - 12.5825;
xCOPRnew2 = xCOPR2 + 12.5825;

% Calculate COP for overall system using Winter et al., 2003 [1]: 
% Equation 1 (same Equation 1 in the Lab 3 manual).

xCOP_f = (xCOPLnew.*(Fz_L./(Fz_L+Rv_R))) + (xCOPRnew.*(Rv_R./(Fz_L+Rv_R)));
yCOP_f = (yCOPL.*(Fz_L./(Fz_L+Rv_R))) + (yCOPR.*(Rv_R./(Fz_L+Rv_R)));

xCOP_s = (xCOPLnew2.*(Fz_L2./(Fz_L2+Rv_R2))) + (xCOPRnew2.*(Rv_R2./(Fz_L2+Rv_R2)));
yCOP_s = (yCOPL2.*(Fz_L2./(Fz_L2+Rv_R2))) + (yCOPR2.*(Rv_R2./(Fz_L2+Rv_R2)));

%%% Filter all data
% Low-pass filter: Butterworth, fc=10Hz
fc=10;
[b, a] = butter(4, fc/(fs/2),'low');
%Filtering 
xCOP_f_filter = filter(b, a, xCOP_f);
yCOP_f_filter = filter(b, a, yCOP_f); 

xCOP_s_filter = filter(b, a, xCOP_s); 
yCOP_s_filter = filter(b, a, yCOP_s);

% Define AP: y-axis; ML: x-axis

AP_f = yCOP_f_filter;
AP_s = yCOP_s_filter;

ML_f = xCOP_f_filter;
ML_s = xCOP_s_filter;

%% Processing EMG data
% Considerer the timing of activation of the Sol. and TA muscles with
% respect to the COP sway in the AP direction.
%Establish gain G
G = 1000;

% V to mV Conversion
sol_f = sol_f.*1000; 
sol_s = sol_s.*1000;

ta_f = ta_f.*1000;
ta_s = ta_s.*1000;

% Define envelope filter parameters, seleting corner frequency of 2.5Hz, butter filter
fc = 2.5;
[b, a] = butter(4, fc/(fs),'low');

% Define time index
t_sol_f = (1:length(sol_f))/fs;
t_sol_s = [1:length(sol_s)]/fs;

t_ta_f = [1:length(ta_f)]/fs;
t_ta_s = [1:length(ta_s)]/fs;

%remove gain
sol_f = sol_f./G; %Soleus Muscle
sol_s = sol_s./G;

ta_f = ta_f./G; % TA Muscle
ta_s = ta_s./G;
 
%normalization
sol_f = sol_f./max(sol_f);    
sol_s = sol_s./max(sol_s); 

ta_f = ta_f./max(ta_f);
ta_s = ta_s./max(ta_s);

%rectification
sol_f = abs(sol_f);
sol_s = abs(sol_s);  

ta_f = abs(ta_f);
ta_s = abs(ta_s);

%signal envelope
sol_f_enveloped = filter(b, a, sol_f);
sol_s_enveloped = filter(b, a, sol_s);
ta_f_enveloped = filter(b, a, ta_f);
ta_s_enveloped = filter(b, a, ta_s);

%Graphing the COP and dual EMG overlay for slow sway
figure;
subplot(2,1,1) 
    hold on
    plot(t_ta_s,ta_s_enveloped)
    plot(t_sol_s,sol_s_enveloped) 
    hold off
    title('EMG: Slow Against Time');
    xlabel('Time (sec)');
    ylabel('Amplitude (mV)');
    legend('TA Muscle','Sol Muscle');
subplot(2,1,2)
    hold on
    plot(t_ta_s,AP_s)
    hold off
    title('COP Net: AP Slow Against Time');
    xlabel('Time (sec)');
    ylabel('Amplitude (mV)');
%Grpahing the COP and dual EMG overlay for fast sway
figure;
subplot(2,1,1) 
    hold on
    plot(t_ta_f,ta_f_enveloped)
    plot(t_sol_f,sol_f_enveloped) 
    hold off
    title('EMG: Fast Against Time');
    xlabel('Time (sec)');
    ylabel('Amplitude (mV)');
    legend('TA Muscle','Sol Muscle');
subplot(2,1,2)
    hold on
    plot(t_ta_f,AP_f)
    hold off
    title('COP Net: AP Fast Against Time');
    xlabel('Time (sec)');
    ylabel('Amplitude (mV)');


