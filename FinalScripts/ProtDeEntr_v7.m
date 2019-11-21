% Script to create normalized power spectra for de-entrainment and
% entrainment conditions. 
% V1: SO 21-3-2018
% script performs the following steps:
%
%   1) read EEG data
%   2) create the bandpass and bandstop filter
%   3) apply the filter
%   4) apply FFT to check how well it worked
%   5) extract multiplication factor (MF)
%   6) normalize within single conditions in the power domain
%   7) equalize the power domain
%   8) append protocol untill length is correct
%   9) save the different protocols as we need it
%   10) save the different protocols as we need it, version 2

clear
dataloc = 'D:\Experiments\VICI\DeEntrainment\ExampleEEG\';
saveloc = 'D:\Experiments\VICI\DeEntrainment\DataStreamerFiles\';
filename = 'S1_EEGbl_pre';
newname = 'S1';
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20180712'));
%addpath('D:\Other\Programs\Matlab\CircStat2012a');

choin = 'P3_Ring';           % channel of interest
IAF = 9.8;               % individual alpha frequency
win = 2;                % window (on one side) extention of IAF
foilim = [1 49];        % limits of pow spectrum that we are interested in
protlength = 60*10;     % in seconds length of final output
ramptimeS = 15;          % in seconds length of ramp time
fsample = 5000;         % sampling rate of data
ramptimeE = 1/fsample;  % in seconds length of ramp time

%% 1) read EEG data
cfg = [];
cfg.dataset = [dataloc filename '.vhdr']; 
cfg.trialdef.triallength = 40; % length of the continueous file you want to read in (from beginning!)
cfg.trialdef.ntrials = 1;
cfgB = ft_definetrial(cfg);

cfg = [];
cfg.dataset = [dataloc filename '.vhdr'];
cfg.trl = cfgB.trl;
cfg.lpfilter = 'yes';
cfg.lpfreq = [200]; % mostly for visualization purposes. Otherwise data looks like crap
cfg.demean = 'yes';
dat = ft_preprocessing(cfg);

cfg =[];
cfg.keeptrials = 'yes';
dat = ft_timelockanalysis(cfg, dat);
[ch chin] = intersect(dat.label, choin); % which channel will be used for analysis

datatoF = squeeze(dat.trial(:,chin,:));
ndatsample = size(datatoF,1);
FreqUse = linspace(0, (floor(ndatsample)/2)./(ndatsample./fsample), floor(ndatsample)/2+1);

%% 2) create the bandpass and bandstop filter
% first a common bandpass filter from 1-49 Hz
normFreq = foilim./(fsample./2);
[bC aC] = butter(2,normFreq, 'pass');

normFreq = [IAF-win IAF+win]./(fsample./2);
[bBP aBP] = butter(2,normFreq, 'pass');
[bBS aBS] = butter(2,normFreq, 'stop');

figure(1)
freqz(bBP,aBP, length(FreqUse));
wl = 0.05;
set(gca, 'xlim', [0 wl]);
subplot(2,1,2)
set(gca,'xlim', [0 wl]);

figure(2)
freqz(bBS,aBS, length(FreqUse));
wl = 0.05;
set(gca, 'xlim', [0 wl]);
subplot(2,1,2)
set(gca,'xlim', [0 wl]);

%% 3) apply the filter
dataOut = filtfilt(bC,aC,datatoF); % forward and reverse filter direction
dataOutBP = filtfilt(bBP,aBP,dataOut);
dataOutBS = filtfilt(bBS,aBS,dataOut);

figure(3)
subplot(2,1,1)
plot([1/fsample:1/fsample:length(dataOut)/fsample], dataOutBP);
subplot(2,1,2)
plot([1/fsample:1/fsample:length(dataOut)/fsample], dataOutBP);
set(gca,'xlim', [0 5]);

figure(4)
subplot(2,1,1)
plot([1/fsample:1/fsample:length(dataOut)/fsample], dataOutBS);
subplot(2,1,2)
plot([1/fsample:1/fsample:length(dataOut)/fsample], dataOutBS);
set(gca,'xlim', [0 5]);

%% 4) apply FFT to check how well it worked
fft_output_tap = fft(dataOut',[], 2);
fft_output_tap_pow = abs(fft_output_tap).^2; % power!
powavg = mean(fft_output_tap_pow(:,1:floor(ndatsample/2)+1),1);

fft_output_tapBP = fft(dataOutBP',[], 2);
fft_output_tap_powBP = abs(fft_output_tapBP).^2; % power!
powavgBP = mean(fft_output_tap_powBP(:,1:floor(ndatsample/2)+1),1);

fft_output_tapBS = fft(dataOutBS',[], 2);
fft_output_tap_powBS = abs(fft_output_tapBS).^2; % power!
powavgBS = mean(fft_output_tap_powBS(:,1:floor(ndatsample/2)+1),1);

FreqUse = linspace(0, (floor(ndatsample)/2)./(ndatsample./fsample), floor(ndatsample)/2+1);
Frin = find(FreqUse > foilim(1) & FreqUse < foilim(2));
Alfr = find(FreqUse > IAF-win & FreqUse < IAF+win);
InFor1F = setdiff(Frin,Alfr);

figure(5)
subplot(3,1,1)
plot(FreqUse(Frin), powavg(Frin));
title('pow EEG data');
ylim = get(gca, 'ylim');
subplot(3,1,2)
plot(FreqUse(Frin), powavgBP(Frin));
title('bandstop');
set(gca, 'ylim', ylim);
subplot(3,1,3)
plot(FreqUse(Frin), powavgBS(Frin));
title('bandpass');
set(gca, 'ylim', ylim);

figure(6)
plot(FreqUse(Frin), powavgBP(Frin)); hold on
plot(FreqUse(Frin), powavgBS(Frin));
legend({'bandpass';'bandstop'});
set(gca, 'xlim', [IAF-win*3 IAF+win*3]);

%% 5) extract multiplication factor (MF)
% complex spectra:
Con = {'All';'Alpha';'DeEntrainment'};
coi = [2 3]; % only these can be used for the MF
C{1} = fft_output_tap;
C{2} = fft_output_tapBP;
C{3} = fft_output_tapBS;

% See which con has maximum MF:
for it = 1:3   
    T{it} = real(ifft(C{it}));
    totpow(it) = sum(abs(C{it}).^2);
    mf(it) = max(abs(T{it}));  % maximum absolute value of time domain 
end
[mfmax in] = max(mf(coi)); % which one is max of the conditions of interest
in = coi(in);

%% 6) normalize within single conditions in the power domain
for it = 1:3  
    C{it} = conv(C{it}, 1./mf(it));
    totpow(it) = sum(abs(C{it}).^2);
    T{it} = real(ifft(C{it}));
end

%% 7) equalize the power domain
[totpowmin in] = min(totpow(coi));
[frin] = coi(in);
for it = 1:3   
    psp = abs(C{it}).^2./(totpow(it)./totpow(frin)); % normalization of power spectrum
    totpowN(it) = sum(psp); % new total power
    C{it} = sqrt(psp).*exp(i.*angle(C{it}));
    T{it} = real(ifft(C{it}));
    DoubleCh(it) = max(abs(T{it}));
    figure(6)
    subplot(3,1,it);
    plot([1/fsample:1/fsample:length(dataOut)/fsample], T{it});
    set(gca, 'xlim', [10 15]);
    title(Con{it});
end

for it = 1:3   
    figure(7)
    subplot(2,3,it);
    plot([1/fsample:1/fsample:length(dataOut)/fsample], T{it});
    title(Con{it});
    subplot(2,3,it+3);
    plot(FreqUse(Frin), abs(C{it}(Frin).^2));
end

%% 8) append untill length is correct
for it = 1:3
   temp = T{it};
   while length(temp) < protlength*fsample
      temp = [temp T{it}];
   end
   T{it} = temp(1:protlength*fsample);
end

%% 9) add triggers and save (version 1)
T{4} = sin(IAF*2*pi*[1/fsample:1/fsample:protlength]);
Con{4} = 'IAF';
newname = 'Version1';

zeroL = 92.5;
TrigL = 0.05;
TrigTimes = [0 100 400];
TrigValue = [2 1 3];
RampSt = 92.5;
RampE = 400;

for it = 1:length(T)
    stimuli = zeros(5,protlength.*fsample);
    for TT = 1:length(TrigTimes)
       stimuli(5,TrigTimes(TT)*fsample+1:TrigTimes(TT)*fsample+fsample*TrigL) = TrigValue(TT); 
    end
    TS = [zeros(1,RampSt*fsample) T{it}];
    % ramp up
    in = RampSt*fsample:RampSt*fsample+ramptimeS*fsample;    
    TS(in) = TS(in).*linspace(0,1, length(in));
    % ramp down
    in = RampE*fsample:RampE*fsample+ramptimeE*fsample;    
    TS(in) = TS(in).*linspace(1,0, length(in));
    TS(in(end)+1:end) = [];
    stimuli = stimuli(:,1:length(TS));
    stimuli(3,:) = TS;
    
    % add to stimuli so that it doesn't start and end at non-zero
    Tt = 0.05;
    stimuli = [zeros(5, Tt*fsample) stimuli zeros(5, Tt*fsample)];    
    save([saveloc newname '_' Con{it}], 'stimuli');
end

%% 10) add triggers and save (version 2)
T{4} = sin(IAF*2*pi*[1/fsample:1/fsample:protlength]);
Con{4} = 'IAF';
newname = 'Version2';

zeroL = 1/fsample;
TrigL = 0.05;
TrigTimes = [7.5 307.5 407.5];
TrigValue = [1 2 3];
RampSt = 1/fsample; 
RampE = 307.5;

for it = 1:length(T)
    stimuli = zeros(5,protlength.*fsample);
    for TT = 1:length(TrigTimes)
       stimuli(5,TrigTimes(TT)*fsample+1:TrigTimes(TT)*fsample+fsample*TrigL) = TrigValue(TT); 
    end
    TS = [zeros(1,RampSt*fsample) T{it}];
    % ramp up
    in = RampSt*fsample:RampSt*fsample+ramptimeS*fsample;    
    TS(in) = TS(in).*linspace(0,1, length(in));
    % ramp down
    in = RampE*fsample:RampE*fsample+ramptimeE*fsample;    
    TS(in) = TS(in).*linspace(1,0, length(in));
    TS(in(end)+1:end) = [];
    maxL = max(length(TS), find(stimuli(5,:) > 0, 1, 'last'));
    stimuli = stimuli(:,1:maxL);
    stimuli(3,1:length(TS)) = TS;
    
    % add to stimuli so that it doesn't start and end at non-zero
    Tt = 0.05;
    stimuli = [zeros(5, Tt*fsample) stimuli zeros(5, Tt*fsample)];    
    save([saveloc newname '_' Con{it}], 'stimuli');
end
plot(stimuli(3,:));
hold on
plot(stimuli(5,:))

