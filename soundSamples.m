
f1 = 500; % low frequency (Hz)
f2 = 2100; % high frequency (Hz)
duration = 1.5; % in seconds
amplitude = 0.2; % normalized between 0 and 1
sampleRate = 20000; % samples per second (Hz)
durationToPlot = 0.01;
writeToFile = 0;
fileExtension = '.wav'; %other options are '.ogg','.flac' and '.mp4'

gatingTime = 0.01; % ramping time at the onset and offset of each sound
gatingSamples = round(gatingTime*sampleRate);
gatingEnveloppe = sqrt(1-cos((0:gatingSamples-1)/(gatingSamples-1)).^2);
samplesToPlot = gatingSamples+(1:round(durationToPlot*sampleRate));

%time vector (in seconds)
t = (1:sampleRate*duration)/sampleRate;

% make low-frequency sinewave
pureTone1 = amplitude * sin(t*f1*2*pi);
%gate onset and offset (to avoid clicks)
pureTone1(1:gatingSamples) = pureTone1(1:gatingSamples) .* gatingEnveloppe;
pureTone1(end-gatingSamples+1:end) = pureTone1(end-gatingSamples+1:end) .* fliplr(gatingEnveloppe);

% make high-frequency sinewave
pureTone2 = amplitude * sin(t*f2*2*pi);
%gate onset and offset (to avoid clicks)
pureTone2(1:gatingSamples) = pureTone2(1:gatingSamples) .* gatingEnveloppe;
pureTone2(end-gatingSamples+1:end) = pureTone2(end-gatingSamples+1:end) .* fliplr(gatingEnveloppe);

% add both tones
bothTones = pureTone1+pureTone2;

%add tones, but shifted in time
timeShift = 0.2; %in seconds
shiftedTones = [pureTone1 zeros(1,timeShift*sampleRate)] + [zeros(1,timeShift*sampleRate) pureTone2];
shiftedT = (1:sampleRate*(duration+timeShift))/sampleRate;
shiftedSamplesToPlot = round((timeShift-durationToPlot/2)*sampleRate):round((timeShift+durationToPlot/2)*sampleRate);

%write files
if writeToFile
  audiowrite(['pureTone1' fileExtension],pureTone1,sampleRate);
  audiowrite(['pureTone2' fileExtension],pureTone2,sampleRate);
  audiowrite(['bothTones' fileExtension],bothTones,sampleRate);
  audiowrite(['shiftedTones' fileExtension],shiftedTones,sampleRate);
end

%play low frequency
sound(pureTone1,sampleRate);
figure;
subplot(4,1,1);
plot(t(samplesToPlot),pureTone1(samplesToPlot));
title(sprintf('Low-frequency tone (%d Hz)',f1));
xlabel('time (s)');
ylim(2*[-amplitude amplitude]);
pause(duration+.5)

%play high frequency
sound(pureTone2,sampleRate);
subplot(4,1,2);
plot(t(samplesToPlot),pureTone2(samplesToPlot));
title(sprintf('High-frequency tone (%d Hz)',f2));
xlabel('time (s)');
ylim(2*[-amplitude amplitude]);
pause(duration+.5)

%play both together
sound(bothTones,sampleRate);
subplot(4,1,3);
plot(t(samplesToPlot),bothTones(samplesToPlot));
title(sprintf('Both tones, simultaneous (%d + %d Hz)',f1,f2));
xlabel('time (s)');
ylim(2*[-amplitude amplitude]);
pause(duration+.5)

%play both, but shifted in time
sound(shiftedTones,sampleRate);
subplot(4,1,4);
plot(shiftedT(shiftedSamplesToPlot),shiftedTones(shiftedSamplesToPlot));
title(sprintf('Both tones, shifted in time (%d + %d Hz)',f1,f2));
xlabel('time (s)');
ylim(2*[-amplitude amplitude]);
