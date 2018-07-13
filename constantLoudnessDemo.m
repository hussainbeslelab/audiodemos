loudnessCorrection = false;
f1=200;
f2=16000;
numberOfFrequencies = 12;

% logarithmically-space frequencies between f1 and f2:
frequencies = 2.^(linspace(log2(f1),log2(f2),numberOfFrequencies)); 
frequencies = [200, 400,  1000, 2000, 5000, 7000, 8000, 9000, 11000, 12000, 14000, 16000];
amplitude = 0.05; % normalized between 0 and 1
sampleRate = 48000; % samples per second (Hz)
duration=2;
pauseDuration=2;
durationToPlot = 0.01;
writeToFile = 1;
fileExtension = '.wav'; %other options are '.ogg','.flac' and '.mp4'

gatingTime = 0.02; % ramping time at the onset and offset of each sound
gatingSamples = round(gatingTime*sampleRate);
gatingEnveloppe = sqrt(1-cos((0:gatingSamples-1)/(gatingSamples-1)).^2);
samplesToPlot = gatingSamples+(1:round(durationToPlot*sampleRate));

% correct for loudness using A-weighting
loudnessConstant = 1;
RA = @(f)12200^2*f.^4./((f.^2+20.6^2).*sqrt((f.^2+107.7^2).*(f.^2+737.9^2)).*(f.^2+12200^2));
A = @(f)2+20*log10(RA(f));

%time vector (in seconds)
t = (1:sampleRate*duration)/sampleRate;

fig=figure;
subplot(2,1,1);
subplot(2,1,2);

soundSequence=[];
for f=frequencies
  % make low-frequency sinewave
  pureTone = amplitude * sin(t*f*2*pi);
  %gate onset and offset (to avoid clicks)
  pureTone(1:gatingSamples) = pureTone(1:gatingSamples) .* gatingEnveloppe;
  pureTone(end-gatingSamples+1:end) = pureTone(end-gatingSamples+1:end) .* fliplr(gatingEnveloppe);

  if loudnessCorrection
    pureTone = pureTone/RA(f)*loudnessConstant;
  end
  %play low frequency
  sound(pureTone,sampleRate);

  if f~=frequencies(1);
    delete(h1);
    delete(h2);
  end
  
  subplot(2,1,1);
  h1 = plot(t(samplesToPlot),pureTone(samplesToPlot));
  hold on;
  title(sprintf('%d Hz tone (Actual amplitude)',round(f)));
  if f==frequencies(1);
    xlabel('time (s)');
    ylim(2*[-amplitude amplitude]);
  end
  
  %plot perceived loudness
  subplot(2,1,2);
  h2=plot(t(samplesToPlot),RA(f)*pureTone(samplesToPlot));
  hold on;
  title(sprintf('%d Hz tone (Perceived amplitude)',round(f)));
  if f==frequencies(1);
    xlabel('time (s)');
    ylim(2*[-amplitude amplitude]);
  end
  
  pause(duration+pauseDuration/2)
  
  subplot(2,1,1);
  delete(h1);
  h1=plot(t(samplesToPlot),zeros(size(samplesToPlot)));
%   ylim(2*[-amplitude amplitude]);
  subplot(2,1,2);
  delete(h2);
  h2=plot(t(samplesToPlot),zeros(size(samplesToPlot)));
%   ylim(2*[-amplitude amplitude]);
  pause(pauseDuration/2)

  soundSequence = [soundSequence pureTone zeros(size(pureTone))];
end


%write files
if writeToFile
  if loudnessCorrection
    audiowrite('constantLoudness.wav',soundSequence,sampleRate);
  else
    audiowrite('varyingLoudness.wav',soundSequence,sampleRate);
  end
end

close(fig);