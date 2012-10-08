% frequencies
Fs = 44100;
Fmin = 100;
Fmax = 10000;
stimdur = 500;
corr_frange = [200 10000];

% plotting limits
dblim = [-100 -20];
freqlim = [0 12];

MIN_DB = -120;


%% load dummy calibration  data
caldata = fake_caldata('FREQS', 100:100:10000);

% assign meaningful values to caldata.mag
x = 2*pi*linspace(0, 1, length(caldata.mag(1, :)));
m = 20 * sin(x);
caldata.mag = [m; m];

% shorter names for freq and mag from caldata
F = caldata.freq;
M = caldata.mag(1, :);

% plot
figure(1)
plot(0.001*F, M, '.-')
xlabel('Frequency (kHz)')
ylabel('dB SPL')
title('System Transfer Function')

%% synthesize test noise from min(F) to max(F)
s = synmononoise_fft(stimdur, Fs, Fmin, Fmax, 1, 0);
s = normalize(s);
s = sin2array(s, 5, Fs);

% s = synmonosine(stimdur, Fs, Fsine, 1, 0);
% plot spectrum of s
[fraw, magraw, phiraw] = daqdbfullfft(s, Fs, length(s));
figure(2)
tvec = 1000 * (0:(length(s)-1)) ./ Fs;
subplot(4, 3, 1)
plot(tvec, s)
title('Test signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(4, 3, 2)
plot(0.001*fraw, magraw);
title('Test Signal Magnitude')
xlabel('freq (kHz)')
ylabel('dB')
ylim(dblim);
xlim(freqlim);

subplot(4, 3, 3)
plot(0.001*fraw, phiraw);
title('Test Signal Phase');
xlabel('freq (kHz)');
ylabel('rad');
xlim(freqlim);;

% play sound
fprintf('Playing raw sound...');
sound(sin2array(0.9*s, 1, Fs), Fs);
pause(1);
fprintf('\n');


%% use compensate signal to compensate with BOOST
[sboost, Sfull, Hboost, Fboost] = compensate_signal(s, F, M, Fs, corr_frange, 'Method', 'BOOST', 'Normalize', 'on', 'Lowcut', 'off');

% plot compensated signal
[fboost, magboost, phiboost] = daqdbfullfft(sboost, Fs, length(sboost));
figure(2)
subplot(4, 3, 4)
plot(tvec, sboost)
title('Boost Comp Signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(4, 3, 5)
plot(0.001*fboost, magboost);
title('Boost Comp Magn.')
xlabel('freq (kHz)')
ylabel('dB')
ylim(dblim);
xlim(freqlim);

subplot(4, 3, 6)
plot(0.001*fboost, phiboost);
title('Boost Comp Phase');
xlabel('freq (kHz)');
ylabel('rad');
xlim(freqlim);;

% play scaled and windowed sound
fprintf('Playing BOOST compensated sound...')
sound(sin2array(0.9*sboost, 1, Fs), Fs);
pause(1);
fprintf('\n');

%% test ATTEN method
% apply correction (ATTEN method)
[satten, Satten, Hatten, foutatten] = compensate_signal(s, F, M, Fs, corr_frange, 'Method', 'ATTEN', 'Normalize', 'on', 'Lowcut', 'off');

figure(2)
% plot compensated signal
[fatten, magatten, phiatten] = daqdbfullfft(satten, Fs, length(satten));
subplot(4, 3, 7)
plot(tvec, satten)
title('Atten Comp Signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(4, 3, 8)
plot(0.001*fatten, magatten);
title('Atten Comp Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim(dblim);
xlim(freqlim);

subplot(4, 3, 9)
plot(0.001*fatten, phiatten);
title('Atten Comp Phase');
xlabel('freq (kHz)');
ylabel('rad');
xlim(freqlim);;

% play sound
fprintf('Playing ATTEN compensated sound...')
sound(sin2array(0.9*satten, 1, Fs), Fs);
fprintf('\n')
pause(1);

%% test COMPRESS method
% apply correction (COMPRESS method)
[scomp, Scomp, Hcomp, foutcomp] = compensate_signal(s, F, M, Fs, corr_frange, 'Method', 'COMPRESS', 'Normalize', 'on', 'Lowcut', 'off');

figure(2)
% plot compensated signal
[fcomp, magcomp, phicomp] = daqdbfullfft(scomp, Fs, length(scomp));
subplot(4, 3, 10)
plot(tvec, scomp)
title('Compress Comp Signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(4, 3, 11)
plot(0.001*fcomp, magcomp);
title('Compress Comp Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim(dblim);
xlim(freqlim);

subplot(4, 3, 12)
plot(0.001*fcomp, phicomp);
title('Compress Comp Phase');
xlabel('freq (kHz)');
ylabel('rad');
xlim(freqlim);;

% play sound
fprintf('Playing COMPRESS compensated sound...')
sound(sin2array(0.9*scomp, 1, Fs), Fs);
pause(1);
fprintf('\n');

