% Script to create normalized power spectra for de-entrainment and
% entrainment conditions. 
% V1: SO 9-3-2018
% script performs the following steps:
%
%   1) read EEG data
%   2) do fft for all points
%   3) extract frequencies of interest
%   4) create 1 over F spectra or the regular EEG 
%   5) interpolate data based on window of interest
%   6) creatrumte phase spec
%   7) create pow spectrum with alpha set at zero
%   8) for comparison, create only alpha pow spectrum
%   9) go to time domain and extract multiplication factor in time comain (MF)
%   10) convolution in the frequency domain for normalization
%   11) still total power is not identical. Final adjustment in freq domain to get the power spectra the same

clear
dataloc = 'D:\Sanne\DeEntrainment\Pilot_TES\Protocol\ExampleEEG\';
saveloc = 'D:\Sanne\DeEntrainment\Pilot_TES\Protocol\DataStreamerFiles\';
filename = 'P001_S4_IAF_pre2';
newname = 'TEST';
addpath('D:\Sanne\DeEntrainment\Pilot_TES\Protocol\fieldtrip-20161231\');
addpath('D:\Sanne\DeEntrainment\Pilot_TES\Protocol\CircStat2012a');

choin = 'Pz';           % channel of interest
IAF = 11.5;             % individual alpha frequency
win = 2;                % window (on one side) extention of IAF
foilim = [1 49];        % limits of pow spectrum that we are interested in
protlength = 60*10;     % in seconds length of final output
ramptime = 10;          % in seconds length of ramp time
Ph = 0;                 % what type of phase estimate (0 = fully random; 1 = only non-estimated random (takes longer); 2 = nearest interpolated)
F1 = 1;                 % Take original EEG or 1 over F (0 = original; 1 = 1 over F);
fsample = 5000;         % sampling rate of data

%% 1) read EEG data
cfg = [];
cfg.dataset = [dataloc filename '.vhdr']; 
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.prestim = 1;
cfg.trialdef.poststim = 1+1/fsample; % we need an odd number so add 1 samplepoint!
cfg.trialdef.eventvalue = {'S  4'; 'S 16'}; % (4,16 is open, 8,32 is closed) > might be different for other data
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
ndatsample = size(datatoF,2);

%% 2) do fft for all points:
tap = hanning(ndatsample)';
tap = tap./norm(tap, 'fro');
fft_output_tap = fft(bsxfun(@times,datatoF,tap),[], 2);
%fft_output_tap = fft_output_tap ./norm(fft_output_tap , 'fro');
fft_output_tap_pow = abs(fft_output_tap).^2; % power!
fft_output_tap_pha = angle(fft_output_tap); % angle
powavg = mean(fft_output_tap_pow(:,1:floor(ndatsample/2)+1),1);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'pow';
%cfg.foilim = [1 49];
%cfg.pad = 0;
cfg.channel = chin;
frft = ft_freqanalysis(cfg, dat);
figure
plot(frft.freq, frft.powspctrm);
set(gca,'xlim',[1 49]);

%% 3) Extract frequencies of interest.
FreqUse = linspace(0, (floor(ndatsample)/2)./(ndatsample./fsample), floor(ndatsample)/2+1);
Frin = find(FreqUse > foilim(1) & FreqUse < foilim(2));
Alfr = find(FreqUse > IAF-win & FreqUse < IAF+win);
InFor1F = setdiff(Frin,Alfr);

figure(1)
plot(FreqUse(Frin), powavg(Frin));
title('pow EEG data');

%% 4) decide what to use for the all 1 over F or the regular EEG? 
if F1 == 1 % do the 1 over F correction
    % fit c/f^a to this data and plot to see how it looks
    clear x
    F = @(x,xdata) x(1)./(xdata.^x(2)); % the function
    LB = [-2 -2];
    UB = [2 2];
    
    for it = 1:30 % to avoid going into local minima start at different x0
        x0 = rand(1,2).*(UB-LB)+LB;
        [xit{it}, resnorm(it)] = lsqcurvefit(F, x0, FreqUse(InFor1F), powavg(InFor1F), LB, UB);
    end;
    [v in] = min(resnorm);
    x = xit{in};
    
    figure(2)
    plot(FreqUse(Frin), powavg(Frin));
    hold on
    plot(FreqUse(InFor1F),F(x,FreqUse(InFor1F)),'-o');
    title('fit');
    
    PowSpNew = zeros(1,ndatsample);
    PowSpNew(Frin) = F(x,FreqUse(Frin));
    PowSpNew(ndatsample-Frin) = F(x,FreqUse(Frin)); % other end of the spectrum
    
    % possible: adjust function to end at zero by subtracting last value:
    % x(3) = F(x,FreqUse(InFor1F(end)));
    % F = (@(x,xdata) x(1)./(xdata.^x(2))-x(3));
else % no do the 1/f power adjustment:    
    PowSpNew = zeros(1,ndatsample);
    PowSpNew(Frin) = mean(fft_output_tap_pow(:,Frin));
    PowSpNew(ndatsample-Frin) = mean(fft_output_tap_pow(:,Frin));
end

%% 5) interpolate data based on window of interest:
ndatsampleNew = protlength*fsample+1; % amount of new samples
% we need an odd number (to fullfill Hermitian symmetry)
if ~mod(ndatsampleNew,2)
    ndatsampleNew = ndatsampleNew-1;
end
FreqUseT = linspace(0, ndatsample./(ndatsample./fsample), ndatsample);
FreqUseNew = linspace(0, ndatsample./(ndatsample./fsample), ndatsampleNew);

% interpolation
PowSpNewInt = interp1(FreqUseT, PowSpNew,FreqUseNew, 'linear');

% to make sure it is really symmetrical
PowSpNewInt(ceil(ndatsampleNew/2)+1:end) = PowSpNewInt(floor(ndatsampleNew/2):-1:1);

figure(3)
plot(FreqUseT, PowSpNew); hold on;
plot(FreqUseNew, PowSpNewInt);

%% 6) create phase spectrum:
% this is to calculate the relative phase:
IAFin = nearest(FreqUse,IAF);
relpha = circ_dist(fft_output_tap_pha, repmat(fft_output_tap_pha(:,IAFin),[1 ndatsample]));
meanrelalph = circ_mean(relpha);

% try to interpolate doesn't make sense:
% PhSpcNewInt = interp1(FreqUseT, meanrelalph,FreqUseNew, 'nearest');

% random phase:
PhSpcNewInt = rand(size(FreqUseNew)).*2*pi-pi;
if Ph == 1
    % original phase for the ones we can estimate
    f = arrayfun(@(x) nearest(FreqUseNew,x),FreqUseT);
    PhSpcNewInt(f) = meanrelalph;
end

% try to interpolate doesn't make sense as phase is random!:
if Ph == 2
    PhSpcNewInt = interp1(FreqUseT, meanrelalph,FreqUseNew, 'nearest');
end

% make symmetric
PhSpcNewInt(ceil(ndatsampleNew/2)+1:end) = PhSpcNewInt(floor(ndatsampleNew/2):-1:1).*-1;
PhSpcNewInt = ifftshift(PhSpcNewInt);

%% 7) create pow spectrum with alpha set at zero:
Alfr = find(FreqUseNew >= IAF-win & FreqUseNew <= IAF+win);
Frin = find(FreqUseNew > foilim(1) & FreqUseNew < foilim(2));
FrinR = setdiff(Frin,Alfr);

PowSpNewIntCut = PowSpNewInt;
PowSpNewIntCut(Alfr) = 0;
PowSpNewIntCut(ndatsampleNew-Alfr) = 0;

figure(4)
plot(FreqUseNew, PowSpNewIntCut);
hold on

% add values: not sure this is needed anymore...:
totpowA = sum(PowSpNewInt(Alfr));
addval = totpowA./length(FrinR);
PowSpNewIntCut(FrinR) = PowSpNewIntCut(FrinR)+addval;
PowSpNewIntCut(ndatsampleNew-FrinR) = PowSpNewIntCut(FrinR)+addval;

plot(FreqUseNew, PowSpNewIntCut);

%% 8) for comparison, create only alpha:
PowSpNewIntA = zeros(size(PowSpNewInt));
PowSpNewIntA(Alfr) = PowSpNewInt(Alfr);
PowSpNewIntA(ndatsampleNew-Alfr) = PowSpNewInt(Alfr);

plot(FreqUseNew, PowSpNewIntA);

%% 9) go to time domain and extract multiplication factor (MF)
% complex spectra:
Con = {'All';'DeEnt';'Alpha'};
coi = [1 2]; % only these can be used for the MF
C{1} = sqrt(PowSpNewInt).*exp(i.*PhSpcNewInt);
C{2} = sqrt(PowSpNewIntCut).*exp(i.*PhSpcNewInt);
C{3} = sqrt(PowSpNewIntA).*exp(i.*PhSpcNewInt);

% See which con has maximum MF:
for it = 1:3   
    T{it} = real(ifft(C{it}));
    mf(it) = max(abs(T{it}));  % maximum absolute value of time domain 
end
[mfmax in] = max(mf(coi)); % which one is max of the conditions of interest
in = coi(in);

%% 10) normalize with this value in the frequency domain, i.e. multi in time is convol in freq: 
for it = 1:3  
    C{it} = conv(C{it}, 1./mfmax);
    totpow(it) = sum(abs(C{it}).^2);
    T{it} = real(ifft(C{it}));
end

%% 11) still totpower is not identical. Final adjustment in Freq domain to get the power spectra the same:
[totpowmin in] = min(totpow(coi));
[frin] = coi(in);
for it = 1:3   
    psp = abs(C{it}).^2./(totpow(it)./totpow(frin)); % normalization of power spectrum
    totpowN(it) = sum(psp); % new total power
    C{it} = sqrt(psp).*exp(i.*PhSpcNewInt);
    T{it} = real(ifft(C{it}));
    DoubleCh(it) = max(abs(T{it}));
    figure(6)
    subplot(3,1,it);
    plot([0:1/fsample:protlength], T{it});
    set(gca, 'xlim', [10 15]);
    title(Con{it});
end

for it = 1:3   
    figure(7)
    subplot(3,1,it);
    plot([0:1/fsample:protlength], T{it});
    title(Con{it});
end

%% now save the different protocols as we need it:
T{4} = sin(IAF*2*pi*[0:1/fsample:protlength]);
Con{4} = 'IAF';

stimuli = zeros(5,protlength.*fsample+1);
stimuli(5,10:floor(fsample.*0.05)) = 255;
stimuli(5,end-floor(fsample.*0.05)-10:end-10) = 255;
for it = 1:length(T)
    stimuli(3,:) = T{it};
    % ramp up
    in = 1:ramptime*fsample;
    stimuli(3,in) = T{it}(in).*linspace(0,1, length(in));
    % ramp down
    in = protlength*fsample-ramptime*fsample:protlength*fsample+1;
    stimuli(3,in) = T{it}(in).*linspace(1,0, length(in));
    save([saveloc newname Con{it}], 'stimuli');
end