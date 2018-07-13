
f_pt = 1000; % frequency of pure tone(Hz)
bandwidth = [inf 20 50 100 200 400 800]; %bandwidth  of noise/notch
bandpass = false; % if set to false, uses notch noise instead of bandlimited noise
pt_duration = .5; % in seconds
noise_duration = 3; % in seconds
pt_amplitude = 0.07; % amplitude of pure tone, normalized between 0 and 1
noiseAmplitude = 0.5; %(rms) amplitude of noise
sampleRate = 20000; % samples per second (Hz)
durationToPlot = 0.1;
maxFrequencyPlot = 2*f_pt;
writeToFile = 0;
fileExtension = '.wav'; %other options are '.ogg','.flac' and '.mp4'
offsetTime = (noise_duration - pt_duration)/2; %when the pure tone is presented within the noise

gatingTime = 0.01; % ramping time at the onset and offset of each sound
gatingSamples = round(gatingTime*sampleRate);
gatingEnveloppe = sqrt(1-cos((0:gatingSamples-1)/(gatingSamples-1)).^2);
samplesToPlot = gatingSamples + offsetTime*sampleRate + (1:round(durationToPlot*sampleRate));

%time vectors (in seconds)
t_pt = (1:sampleRate*pt_duration)/sampleRate;
nSamples_pt = length(t_pt);
t_noise = (1:sampleRate*noise_duration)/sampleRate;
nSamples_noise = length(t_noise);

% make pure tone
pureTone = pt_amplitude * sin(t_pt*f_pt*2*pi);
%gate onset and offset (to avoid clicks)
pureTone(1:gatingSamples) = pureTone(1:gatingSamples) .* gatingEnveloppe;
pureTone(end-gatingSamples+1:end) = pureTone(end-gatingSamples+1:end) .* fliplr(gatingEnveloppe);

figure('units','normalized','outerposition',[0 0 1 1]);
nBandwidths=length(bandwidth);
for b = 1:nBandwidths
  % make noise
  noise = rand(size(t_noise))-0.5;
  noise = noise/rms(noise)*noiseAmplitude; %set rms level 
  %gate onset and offset (to avoid clicks)
  noise(1:gatingSamples) = noise(1:gatingSamples) .* gatingEnveloppe;
  noise(end-gatingSamples+1:end) = noise(end-gatingSamples+1:end) .* fliplr(gatingEnveloppe);
  % filter noise
  noise_spectrum = fft(noise);
  noise_phase = angle(noise_spectrum);
  noise_amplitude = abs(noise_spectrum);
  % select frequencies
  bp_low = f_pt - bandwidth(b)/2; % low frequency of bandpass/notch noise
  bp_high = f_pt + bandwidth(b)/2; % high frequency of bandpass/notch noise
  frequencies_noise = sampleRate*(0:nSamples_noise/2-1)/nSamples_noise;
  if  bandpass
    filter = (frequencies_noise < bp_high) & (frequencies_noise > bp_low);
  else
    filter = (frequencies_noise > bp_high) | (frequencies_noise < bp_low);
  end
  %duplicate filter for negative frequencies
  filter = [filter fliplr(filter)];
  %apply filter
  noise_amplitude = noise_amplitude .* filter;
  noise_spectrum = noise_amplitude .* exp(1j * noise_phase); 
  % inverse FFT
  noise = real(ifft(noise_spectrum));
  % noise = noise/rms(noise)*noiseAmplitude; %set rms level 

  % add both sound
  bothSounds = noise;
  bothSounds(offsetTime*sampleRate+(1:pt_duration*sampleRate)) = pureTone+noise(offsetTime*sampleRate+(1:pt_duration*sampleRate));
  both_spectrum = fft(bothSounds);

  % %add tones, but shifted in time
  % timeShift = 0.2; %in seconds
  % shiftedTones = [pureTone zeros(1,timeShift*sampleRate)] + [zeros(1,timeShift*sampleRate) noise];
  % shiftedT = (1:sampleRate*(duration+timeShift))/sampleRate;
  % shiftedSamplesToPlot = round((timeShift-durationToPlot/2)*sampleRate):round((timeShift+durationToPlot/2)*sampleRate);

  %write files
  if writeToFile
    audiowrite(['bothSounds' fileExtension],bothSounds,sampleRate); %need to change the file name
  end

  %play both together
  subplot(nBandwidths,2,(b-1)*2+1);
  plot(t_noise(samplesToPlot),bothSounds(samplesToPlot));
  if b==1
    title('Pure tone embedded in noise (time domain)');
  end
  xlabel('time (s)');
  ylim(2*[-noiseAmplitude noiseAmplitude]);
  subplot(nBandwidths,2,(b-1)*2+2)
  plot(frequencies_noise,abs(both_spectrum(1:end/2)));
  if b==1
    title('Spectrum');
  end
  xlabel('Frequency (Hz)');
  xlim([0 maxFrequencyPlot]);
  sound(bothSounds,sampleRate);
  pause(noise_duration+.5)

end
