function SMS()

%==============================================================
% SMS-Matlab like emulation
%%==============================================================
clear all
close all
%==== USER DATA =====
DAFx_in = audioread('AfluteNaFr-C5.wav'); % wave file
SR = 22050; % sampling rate
w1Length = 2048; % analysis window size
n1 =512 ; % analysis window hop size
nPeaks = 100; % number of peaks detected
nSines = 100; % number of sinuosoids to track (and synthetise)
minSpacePeaks = 2; % minimum space (bins) between two picked peaks
zp = 2; % zero-padding coefficient
rgain = 1.; % gain for the residual component
MaxFreq = 11000; % maximum frequency, in Hertz, for plottings
MinMag = -100; % minimum magnitude, in dB, for plottings
%--- figure data
fig1 = 'yes'; % if uncommented, will plot the Blackman-Harris
% window
fig2 = 'yes'; % if uncommented, will plot the peaks detection
% and tracking in one frame
fig3 = 'yes'; % if uncommented, will plot the peak trackings
% real-time
fig4 = 'yes'; % if uncommented, will plot the original and
% the transformed FFT in one frame
fig5 = 'yes'; % if uncommented, will plot the peak trackings
% only at the end of the process
fig6 = 'yes'; % if uncommented, will plot the original signal,
% its sine and residual part, and the transformed signal
%=== Definition of the Windows ===
%--- definition of the analysis window
fConst=2*pi/(w1Length+1-1);
w1=[1:w1Length]';
w1=.35875 -.48829*cos(fConst*w1)+.14128*cos(fConst*2*w1) ...
    -.01168*cos(fConst*3*w1);
w1=w1/sum(w1)*2;
N=w1Length*zp; % new size of the window
%--- synthesis window
w2=w1;
n2=n1;
%--- triangular window
wt2=triang(n2*2+1); % triangular window
%--- main lobe table of bh92
[bh92SINE2SINE,bh92SINE2SINEsize]=bh92SINE2SINEgeneration;
%--- data for the loops
frametime = n1/SR;
pin = 0;
pout = 0;
TuneLength=length(DAFx_in);
pend=TuneLength-w1Length;
%=== Definition of the data arrays ===
DAFx_in = [zeros(w1Length/2-n1-1,1); DAFx_in];
DAFx_outsine = zeros(TuneLength,1);
DAFx_outres = zeros(TuneLength,1);
%--- arrays for the partial tracking
iloc = zeros(nSines,1);
ival = zeros(nSines,1);
iphase = zeros(nSines,1);
previousiloc = zeros(nSines,1);
previousival = zeros(nSines,1);
maxSines = 400; % maximum voices for harmonizer
syniloc = zeros(maxSines,1);
synival = zeros(maxSines,1);
previoussyniloc = zeros(maxSines,1);
previousiphase = zeros(maxSines,1);
currentiphase = zeros(maxSines,1);
%--- arrays for the sinus' frequencies and amplitudes
SineFreq = zeros(nSines,ceil(TuneLength/n2));
SineAmp = zeros(nSines,ceil(TuneLength/n2));
pitch = zeros(1,1+ceil(pend/n1));
pitcherr = zeros(1,1+ceil(pend/n1));
%--- creating figures ---
if(exist('fig1'))
    h = figure(1); set(h,'position', [10, 45, 200, 200]);
end
if(exist('fig2'))
    h = figure(2); set(h,'position', [10, 320, 450, 350]);
    axisFig2 = [0 MaxFreq MinMag 0]; zoom on;
end
if(exist('fig3'))
    h = figure(3); set(h,'position', [220, 45, 550, 200]);
    axisFig3 = [1 1+ceil(pend/n1) 0 MaxFreq]; zoom on;
end
if(exist('fig4'))
    h = figure(4); set(h,'position', [470, 320, 450, 350]);
    axisFig4 = [0 MaxFreq MinMag 0]; zoom on;
end
if(exist('fig5'))
    h = figure(5); set(h,'position', [220, 45, 550, 200]);
    axisFig5 = [1 1+ceil(pend/n1) 0 MaxFreq]; zoom on;
end
%--- plot the Blackman-Harris window
if(exist('fig1'))
    figure(1)
    plot(20*log10(abs(fftshift(fft(bh92SINE2SINE)/bh92SINE2SINEsize))))
    title('Blackman-Harris window');xlabel('Samples');
    ylabel('Amplitude')
end
tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
disp('analazing frame ...');
while pin<pend
    %--- windowing
    grain = DAFx_in(pin+1:pin+w1Length).*w1(1:w1Length);
    %--- zero padding
    padgrain = zeros(N,1);
    padgrain(1:w1Length/2) = grain(w1Length/2+1:w1Length); %zero phase
    padgrain(N-w1Length/2+1:N) = grain(1:w1Length/2);
    %--- fft computation
    f = fft(padgrain);
    r = abs(f);
    phi = angle(f);
    ft = r.*exp(j*phi);
    %===== Analysis =====
    %--- peak detection (and their plottings)
    [ftloc, ftval]=PickPeaks(r(1:N/2),nPeaks,minSpacePeaks);
    %--- calculate interpolated values (peak position,phase,amplitude)
    [iftloc, iftphase, iftval] = ...
    interpolatedValues (r,phi,N,zp,ftloc,ftval);
    %--- pitch detection
    [pitchvalue,pitcherror,isHarm] = ...
    pitchDetection (r,N,SR,nPeaks,iftloc,iftval);
    pitch(1+pin/n1) = pitchvalue*isHarm;
    pitcherr(1+pin/n1) = pitcherror;
    %--- peaks tracking
    if (pin==0) %--- for the first frame
        nNewPeaks = nSines;
    else %--- creating new born tracks
        for i=1:nSines
            if (previousiloc(i)==0)
            [previousiloc(i), previousival(i)] = CreateNewTrack ...
                (iftloc, iftval, previousiloc, previousival, nSines, MinMag);
            nNewPeaks = nNewPeaks - 1;
            end
        end
        %--- simple Peak tracker
        [iloc,ival,iphase,previousiloc,previousival,distminindex] = ...
        peakTrackSimple(nSines,nPeaks,N,SR,pitchvalue,iftloc, ...
        iftval,iftphase,isHarm,previousiloc,previousival);
    end
    %--- savings
    previousival = ival;
    previousiloc = iloc;
    SineFreq(:,1+pin/n1)=max((iloc-1)/N*SR,0.);
    % frequency of the partials
    SineAmp(:,1+pin/n1)=max(ival, MinMag);
    %==========> amplitudes of the partials
    syniloc(1:nSines) = max(1,iloc);
    synival(1:nSines) = ival;
    if(exist('fig3')) % plot: the trackings of partials
        figure(3); clf; hold on
        PlotTracking(SineFreq(:,1:1+pin/n1), pitch(1:1+pin/n1));
        xlabel('Frame number');ylabel('Frequency (Hz)');
        axis(axisFig3);title('Peak tracking'); drawnow
    end
    %--- residual computation
    resfft = ft;
    if(isHarm==1)
        %alex
        %resfft=resfft-sinefillSpectrum(iloc,ival,iphase,nSines,...
        resfft=resfft-sinefillspectrum(iloc,ival,iphase,nSines,...
            w1Length, zp, bh92SINE2SINE, bh92SINE2SINEsize);
    end
    %--- figures
    if(exist('fig2'))
        figure(2); clf; hold on
        % plot: FFT of the windowed signal (Hz,dB)
        plot((1:N/2)/N*SR, 20*log10(r(1:N/2)));
        for l=1:nPeaks % plot: the peaks detected
            plot([ftloc(l)-1 ftloc(l)-1]/N*SR, ...
            [20*log10(ftval(l)),MinMag-1],'r:x');
        end
        for l=1:nSines % plot: sines tracked and the residual part
            plot([iloc(l)-1, iloc(l)-1]/N*SR, [ival(l), MinMag-1],'k')
        end
        plot((1:N/2)/N*SR, 20*log10(abs(resfft(1:N/2))),'g');
        if(isHarm) % plot: true pitch of each harmonic
            for l=1:nSines
                plot([pitchvalue*l, pitchvalue*l],[1, MinMag-1],'y:')
            end
        end
        xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');axis(axisFig2);
        title('Peak detection and tracking for one frame'); drawnow
    end
    nSynSines = nSines;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===== Transformations =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===== Synthesis =====
%--- phase computation
    if (pin > 0)
        for i=1:nSynSines
            if (syniloc(i)~=0)
                ifreq = (previoussyniloc(distminindex(i))+ syniloc(i))/2;
                % average bin
                freq = (ifreq-1)/N*SR; % freq in Hz (if loc=1 --> freq=0)
                currentiphase(i)=unwrap2pi(previousiphase(distminindex(i))+...
                    2*pi*freq*frametime);
            end
        end
    end
    previoussynival = synival;
    previoussyniloc = syniloc;
    previousiphase = currentiphase;
    %--- compute sine spectrum
    padsynthft=sinefillspectrum(syniloc, synival, currentiphase,nSynSines,w1Length, zp, bh92SINE2SINE, bh92SINE2SINEsize);
    if (isHarm==0)
        padsynthft = zeros(size(padsynthft));
    end
    %--- residual computation
    respadgrain=real(ifft(resfft));
    resgrain=[respadgrain(N-w1Length/2+1:N); ...
    respadgrain(1:w1Length/2)]./w2(1:w1Length);
    ressynthgrain=wt2(1:n2*2).*resgrain(w1Length/2-n2:w1Length/2+n2-1);
    DAFx_outres(pout+1:pout+n2*2)=DAFx_outres(pout+1:pout+n2*2)+ ...
        ressynthgrain;
    %--- sinusoidal computation
    sinpadgrain=real(ifft(padsynthft));
    singrain=[sinpadgrain(N-w1Length/2+1:N); ...
    sinpadgrain(1:w1Length/2)]./w2(1:w1Length);
    sinsynthgrain=wt2(1:n2*2).*singrain(w1Length/2-n2:w1Length/2+n2-1);
    DAFx_outsine(pout+1:pout+n2*2)=DAFx_outsine(pout+1:pout+n2*2)+ ...
    sinsynthgrain;
    %--- figure with original signal and transformed signal FFT
    synthr = abs(fft(respadgrain + sinpadgrain));
    if(exist('fig4'))
        figure(4); clf; hold on
        plot((1:N/2)/N*SR, 20*log10(r(1:N/2)),'b:'); axis(axisFig4);
        plot((1:N/2)/N*SR, 20*log10(synthr(1:N/2)),'r');
        figure(4);
        xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');axis(axisFig4);
        title('FFT of the original (blue) and the')
        title('transformed (red) signals');
        drawnow
    end
    %--- increment loop indexes
    pin = pin + n1;
    pout = pout + n2;
    disp(pin/n1);
end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc
%===== write output sounds =====
DAFx_in = DAFx_in(w1Length/2-n1:length(DAFx_in));
% remove the zeros added for the process
DAFx_outresynth = DAFx_outsine(1:TuneLength)+ ...
rgain*DAFx_outres(1:TuneLength);
mm = max(abs(DAFx_outresynth));
audiowrite('DAFx_out.wav',DAFx_outresynth/mm, SR);
audiowrite('DAFx_outsine.wav',DAFx_outsine/mm, SR);
audiowrite('DAFx_outres.wav',DAFx_outres/mm, SR);
if(exist('fig3')==0 & exist('fig5'))% plot: trackings of partials
    % only at the end of the process
    figure(5); clf; hold on
    PlotTracking(SineFreq(:,1:1+pend/n1), pitch(1:1+pend/n1));
    xlabel('Frame number'); ylabel('Frequency (Hz)'); axis(axisFig5);
    title('Peak tracking'); drawnow
end
if(exist('fig6')) % plot the input signal, its sinus
    % and its residual part, and the transformed signal
    figure(6)
    subplot(4,1,1); plot(DAFx_in); xlabel('input signal');
    subplot(4,1,2); plot(DAFx_outsine);xlabel('sinus part');
    subplot(4,1,3); plot(DAFx_outres);xlabel('residual part');
    subplot(4,1,4); plot(DAFx_outresynth);
    xlabel('resynthetized signal');
end

%===================================================
%===================================================
function[bh92SINE2SINE,bh92SINE2SINEsize]=bh92SINE2SINEgeneration;
%function[bh92SINE2SINE,bh92SINE2SINEsize]=bh92SINE2SINEgeneration;
%
% ==> generation of the Blackman-Harris window
% output data:
% bh92SINE2SINEsize: size of the window
% bh92SINE2SINE: (sampled) window
bh92SINE2SINEsize = 4096;
bh92SINE2SINE = zeros(bh92SINE2SINEsize,1);
bh92N = 512;
bh92const = [.35875, .48829, .14128, .01168];
bh92Theta = -4*2*pi/bh92N;
bh92ThetaIncr = 8*2*pi/bh92N/bh92SINE2SINEsize;
for i=1:bh92SINE2SINEsize
    for m=0:3
        bh92SINE2SINE(i)=bh92SINE2SINE(i)-bh92const(m+1)/2*...
            (sine2sine(bh92Theta-m*2*pi/bh92N,bh92N)+...
        sine2sine(bh92Theta+m*2*pi/bh92N,bh92N));
    end;
    bh92Theta = bh92Theta + bh92ThetaIncr;
end;
bh92SINE2SINE = bh92SINE2SINE/bh92SINE2SINE(bh92SINE2SINEsize/2+1);

%====================================================================
%====================================================================
function x = sine2sine( x , N )
% sine2sine function !!!
x = sin((N/2)*x) / sin(x/2);

%====================================================================
%====================================================================
function [loc, val] = PickPeaks(spectrum, nPeaks, minspace)
%function [loc, val] = pickpeaks(spectrum, nPeaks, minspace)
%
%==> peaking the nPicks highest peaks in the given spectrum
% from the greater to the lowest
% data:
% loc: bin number of peaks (if loc(i)==0, no peak detected)
% val: amplitude of the given spectrum
% spectrum: spectrum (abs(fft(signal))
% nPicks: number of peaks to pick
% minspace: minimum of space between two peaks
[r, c] = size(spectrum);
rmin = min(spectrum) - 1;
% ---find a peak, zero out the data around the peak, and repeat
val = ones(nPeaks,c)*-100;
loc = zeros(nPeaks,c);
for k=1:c %--- find all local peaks
    difference = diff([rmin; spectrum(:,k); rmin]); % derivate
    iloc = find(difference(1:r)>= 0 & difference(2:r+1) <= 0);
    % peak locations
    ival = spectrum(iloc,k); % peak values
for p=1:nPeaks
    [val(p,k),l] = max(ival); % find current maximum
    loc(p,k) = iloc(l); % save value and location
    ind = find(abs(iloc(l)-iloc) > minspace);
    % find peaks which are far away
if (isempty(ind))
    break % no more local peaks to pick
end
    ival = ival(ind); % shrink peak value and location array
    iloc = iloc(ind);
end
end

%====================================================================
%====================================================================
function [iftloc, iftphase, iftval] = interpolatedValues ...
(r, phi, N, zp, ftloc, ftval)
%function [iftloc, iftphase, iftval] = interpolatedValues ...
% (r, phi, N, zp, ftloc, ftval)
%
%==> calculatus of the interpolated values of location (bin),
% phase and magnitude by cubic interpolation
% data:
% iftloc: interpolated location (bin)
% iftval: interpolated magnitude
% iftphase: interpolated phase
% ftloc: peak locations (bin)
% ftval: peak magnitudes
% r: modulus of the FFT
% phi: phase of the FFT
% N: size of the FFT
% zp: zero-padding multiplicative coefficient
%--- calculate interpolated peak position in bins (iftloc) ------
leftftval = r((ftloc-1).*((ftloc-1)>0)+((ftloc-1)<=0).*1);
rightftval= r((ftloc+1).*((ftloc+1)<N/2)+((ftloc+1)>=N/2).*(N/2));
leftftval = 20*log10(leftftval);
rightftval= 20*log10(rightftval);
ftval = 20*log10(ftval);
iftloc = ftloc + .5*(leftftval - rightftval) ./ ...
(leftftval - 2*ftval + rightftval);
%--- interpolated ftloc -----------------------------------------
iftloc = (iftloc>=1).*iftloc + (iftloc<1).*1;
iftloc = (iftloc>N/2+1).*(zp/2+1) + (iftloc<=N/2+1).*iftloc;
%--- calculate interpolated phase (iphase) ----------------------
leftftphase = phi(floor(iftloc));
rightftphase= phi(floor(iftloc)+1);
intpfactor = iftloc-ftloc;
intpfactor = (intpfactor>0).*intpfactor ...
+(intpfactor<0).*(1+intpfactor);
diffphase = unwrap2pi(rightftphase-leftftphase);
iftphase = leftftphase+intpfactor.*diffphase;
%--- calculate interpolate amplitude (iftval) -------------------
iftval = ftval-.25*(leftftval-rightftval).*(iftloc-ftloc);

%====================================================================
%====================================================================
function argunwrap = unwrap2pi (arg)
% function argunwrap = unwrap2pi (arg)
%
%==> unwrapping of the phase, in [-pi, pi]
% arg: phase to unwrap
arg = arg - floor(arg/2/pi)*2*pi;
argunwrap = arg - (arg>=pi)*2*pi;

%====================================================================
%====================================================================
function[pitchvalue,pitcherror,isHarm]=pitchDetection(r,N, ...
SR,nPeaks,iftloc,iftval)
% function[pitchvalue,pitcherror,isHarm]= ...
% pitchDetection(r,N,SR,nPeaks,iftloc,iftval)
%
%==> pitch detection function, using the Two-Way Mismatch
% algorithm (see TWM.m)
%
% data:
% r: FFT magnitude
% N: size of the FFT
% SR: sampling rate
% nPeaks: number of peaks tracked
% iftloc, iftval: location (bin) and magnitude of the peak
%--- harmonicity evaluation of the signal
highenergy = sum(r(round(5000/SR*N):N/2)); % 5000 Hz to SR/2 Hz
lowenergy = sum(r(round(50/SR*N):round(2000/SR*N)));
% 50 Hz to 2000 Hz
isHarm = max(0,(highenergy/lowenergy < 0.6));
if (isHarm==1) %-- 2-way mismatch pitch estimation when harmonic
    npitchpeaks = min(50,nPeaks);
    [pitchvalue,pitcherror] = ...
    TWM(iftloc(1:npitchpeaks),iftval(1:npitchpeaks),N,SR);
else
    pitchvalue = 0;
    pitcherror = 0;
end;
%--- in case of two much pitch error,
% signal supposed to be inhamonic
isHarm = min (isHarm,(pitcherror<=1.5));

%====================================================================
%====================================================================
function [pitch, pitcherror] = TWM (iloc, ival, N, SR)
%function [pitch, pitcherror] = TWM (iloc, ival, N, SR)
%
% => Two-way mismatch error pitch detection
% using Bauchamp & Maher algorithm
%
% data:
% iloc: location (bin) of the peaks
% ival: magnitudes of the peaks
% N: number of peaks
% SR: sampling rate
ifreq = (iloc-1)/N*SR; % frequency in Hertz
%--- avoid zero frequency peak
[zvalue,zindex] = min(ifreq);
if (zvalue==0)
    ifreq(zindex) = 1;
    ival(zindex) = -100;
end
ival2 = ival;
[MaxMag,MaxLoc1] = max(ival2);
ival2(MaxLoc1) = -100;
[MaxMag2,MaxLoc2]= max(ival2);
ival2(MaxLoc2) = -100;
[MaxMag3,MaxLoc3]= max(ival2);
%--- pitch candidates
nCand = 10; % number of candidates
pitchc = zeros(1,3*nCand);
pitchc(1:nCand)=(ifreq(MaxLoc1)*ones(1,nCand))./((nCand ...
+1-[1:nCand]));
pitchc(nCand+1:nCand*2)=(ifreq(MaxLoc2)*ones(1,nCand))./ ((nCand ...
+1-[1:nCand]));
pitchc(nCand*2+1:nCand*3)=(ifreq(MaxLoc3)*ones(1,nCand))./((nCand ...
+1-[1:nCand]));
%pitchc=100:300;
harmonic = pitchc;
%--- predicted to measured mismatch error
ErrorPM = zeros(fliplr(size(harmonic)));
MaxNPM = min(10,length(iloc));
for i=1:MaxNPM
    difmatrixPM = harmonic' * ones(size(ifreq))';
    difmatrixPM = abs(difmatrixPM ...
    -ones(fliplr(size(harmonic)))*ifreq');
    [FreqDistance,peakloc] = min(difmatrixPM,[],2);
    Ponddif = FreqDistance .* (harmonic'.^(-0.5));
    PeakMag = ival(peakloc);
    MagFactor = max(0, MaxMag - PeakMag + 20);
    MagFactor = max(0, 1.0 - MagFactor/75.0);
    ErrorPM = ErrorPM ...
        +(Ponddif+MagFactor.*(1.4*Ponddif-0.5));
    harmonic = harmonic+pitchc;
end
%--- measured to predicted mismatch error
ErrorMP = zeros (fliplr(size(harmonic)));
MaxNMP = min(10,length(ifreq));
for i=1:length(pitchc)
    nharm = round(ifreq(1:MaxNMP)/pitchc(i));
    nharm = (nharm>=1).*nharm + (nharm<1);
    FreqDistance = abs(ifreq(1:MaxNMP) - nharm*pitchc(i));
    Ponddif = FreqDistance.* (ifreq(1:MaxNMP).^(-0.5));
    PeakMag = ival(1:MaxNMP);
    MagFactor = max(0,MaxMag - PeakMag + 20);
    MagFactor = max(0,1.0 - MagFactor/75.0);
    ErrorMP(i) = sum(MagFactor.*(Ponddif ...
        +MagFactor.*(1.4*Ponddif-0.5)));
end
%--- total error
Error = (ErrorPM/MaxNPM) + (0.3*ErrorMP/MaxNMP);
[pitcherror, pitchindex] = min(Error);
pitch = pitchc(pitchindex);

%====================================================================
%====================================================================
function[iloc,ival,iphase,previousiloc,previousival, ...
distminindex]=peakTrackSimple(nSines,nPeaks,N, ...
SR,pitchvalue,iftloc,iftval,iftphase,isHarm, ...
previousiloc,previousival);
% function[iloc,ival,iphase,previousiloc,previousival, ...
% distminindex]=peakTrackSimple(nSines,nPeaks,N, ...
% SR,pitchvalue,iftloc,iftval,iftphase,isHarm, ...
% previousiloc,previousival);
%
%==> simplest partial tracking
% data:
% iloc,ival,iphase: location (bin), magnitude
% and phase of peaks (current frame)
% previousiloc,previousival,previousiphase: idem for
% previous frame
% iftloc, iftval, iftphase: idem of all of the peaks in the FT
% distminindex: indexes of the minimum distance
% between iloc and iftloc
% nPeaks: number of peaks detected
% nSines: number of peaks tracked
% N: size of the FFT
% SR: sampling rate
% pitchvalue: estimated pitch value
% isHarm: indicator of harmonicity
tmpharm = pitchvalue; %--- temporary harmonic
iloc = zeros(nSines,1);
MindB = -100;
ival = zeros(nSines,1) + MindB;
iphase = zeros(nSines,1);
distminindex = zeros(nSines,1);
Delta = 0.01;
for i=1:nSines %--- for each sinus detected
    if (isHarm==1) %--- for a harmonic sound
        [closestpeakmag,closestpeakindex]=min(abs((iftloc-1)/N*SR-tmpharm));
        tmpharm = tmpharm + pitchvalue;
    else %--- for an inharmonic sound
    [closestpeakmag,closestpeakindex]=min(abs(iftloc-previousiloc(i)));
    end
    iloc(i) = iftloc(closestpeakindex); %--- bin of the closest
    ival(i) = iftval(closestpeakindex);
    iphase(i) = iftphase(closestpeakindex);
    dist = abs(previousiloc-iloc(i));
    [distminval, distminindex(i)] = min(dist);
end

%====================================================================
%====================================================================
function[newiloc,newival]=CreateNewTrack(iftloc,iftval, ...
previousiloc,previousival,nSines,MinMag);
% function[newiloc,newival]=CreateNewTrack(iftloc,iftval, ...
% previousiloc,previousival,nSines,MinMag);
%
%==> creation of a new track by looking for a new significant
% peak not already tracked
% data: iftlov, iftval: bin number & magnitude of peaks detected
% previousiloc,
% previousival: idem for previous peaks detected
% nSines: number of sines
% MinMag: minimum magnitude (-100 dB) for
% 0 amplitude
%--- removing peaks already tracked
for i=1:nSines
    [min1, ind] = min(abs(iftval - previousival(i)));
    iftval(ind) = MinMag;
end
%--- keeping the maximum
[newival, ind] = max(iftval);
newiloc = iftloc(ind);

%====================================================================
%====================================================================
function PlotTracking(SineFreq, pitch)
%function PlotTracking(SineFreq, pitch)
%
%==> plot the partial tracking
% data:
% SineFreq: frequencies of the tracks
% pitch: frequency of the pitch
[nSines, nFrames] = size(SineFreq);
for n=1:nSines
    f=1;
    while (f<=nFrames)
        while (f<=nFrames & SineFreq(n,f)==0)
            f = f+1;
        end
        iStart = min(f,nFrames);
        while (f<=nFrames & SineFreq(n,f)>0)
            f = f+1;
        end
        iEnd = min(max(1,f-1),nFrames);
        if (iEnd > iStart)
            line((iStart:iEnd), SineFreq(n,iStart:iEnd));
        end
    end
end
h = line((1:nFrames), pitch(1:nFrames));
set(h,'linewidth', 2, 'Color', 'black');

%====================================================================
%====================================================================
function padsynthft=sinefillspectrum(iloc,ival,iphase,nSines, ...
w1Length, zp, bh92SINE2SINE, bh92SINE2SINEsize)
%function padsynthf =sinefillspectrum(iloc,ival,iphase,nSines, ...
% w1Length, zp, bh92SINE2SINE, bh92SINE2SINEsize)
%
%=> compute the spectrum of all the sines in the frequency
% domain, in order to remove it from the signal
% data:
% padsynth:
% iloc, ival, iphase: location (bin), magnitude value (dB)
% and phase of a peak
% nSines: number of sines (=length of ival and iloc)
% w1Length: size of the analysis window
% zp: zero-padding multiplicative coefficient
% bh92SINE2SINE: Blackman-Harris window
% bh92SINE2SINEsize: Blackman-Harris window size
peakmag=10.^(ival/20); % magnitude (in [0;1])
halflobe=8*zp/2-1; % bin number of the half lobe
firstbin=floor(iloc)-halflobe; % first bin for filling positive
% frequencies
firstbin2=floor(w1Length*zp-iloc+2)-halflobe;
% idem for negative frequencies
binremainder=iloc-floor(iloc);
sinphase=sin(iphase);
cosphase=cos(iphase);
findex=1-binremainder;
bh92SINE2SINEindexes =zeros(8*zp,1);
sinepadsynthft=zeros(w1Length*zp+halflobe+halflobe+1,1);
padsynthft =zeros(w1Length*zp,1);
%--- computation of the complex value
for i=1:nSines %--- for each sine
    if (iloc(i)~=0) %--- JUST WORK WITH NON ZEROS VALUES OF iloc !!!
        % -> tracked sines
        beginindex = floor(0.5 + findex(i)*512/zp)+1;
        % alex
        bh92SINE2SINEindexes = [beginindex:512/zp:beginindex+512/zp*(8*zp-1)]';
        %bh92SINE2SINEindexes=[beginindex:512/zp:beginindex ...
        %+512/zp*(8*zp-1)]';
        if (bh92SINE2SINEindexes(8*zp)>bh92SINE2SINEsize)
            bh92SINE2SINEindexes(8*zp)=bh92SINE2SINEsize;
        end
        magsin=bh92SINE2SINE(bh92SINE2SINEindexes) ...
            .*sinphase(i)*peakmag(i);
        magcos=bh92SINE2SINE(bh92SINE2SINEindexes) ...
            .*cosphase(i)*peakmag(i);
        %--- fill positive frequency
        sinepadsynthft(firstbin(i)+halflobe:firstbin(i) ...
            +halflobe+8*zp-1)= ...
        sinepadsynthft(firstbin(i)+halflobe:firstbin(i)+ ...
        halflobe+8*zp-1)+(magcos+j*magsin);
        %--- fill negative frequency
        if (firstbin2(i)+halflobe <= w1Length*zp)
            sinepadsynthft(firstbin2(i)+halflobe:firstbin2(i) ...
                +halflobe+8*zp-1)= ...
            sinepadsynthft(firstbin2(i)+halflobe:firstbin2(i)+ ...
            halflobe+8*zp-1)+(magcos-j*magsin);
        end
    end
end
%--- fill padsynthft
padsynthft=padsynthft+sinepadsynthft(halflobe+1:halflobe+1 ...
+w1Length*zp-1);
padsynthft(1:halflobe) = padsynthft(1:halflobe) + ...
sinepadsynthft(w1Length*zp+1:w1Length*zp+halflobe);
padsynthft(w1Length*zp-halflobe+1:w1Length*zp) = ...
padsynthft(w1Length*zp-halflobe+1:w1Length*zp) ...
+ sinepadsynthft(1:halflobe);

%====================================================================
%====================================================================
function w = triang(n)
% TRIANG Triangular window.
if rem(n,2)
    % It's an odd length sequence
    w = 2*(1:(n+1)/2)/(n+1);
    w = [w w((n-1)/2:-1:1)]';
else
    % It's even
    w = (2*(1:(n+1)/2)-1)/n;
    w = [w w(n/2:-1:1)]';
end

%====================================================================
%====================================================================
function [syniloc synival] = transform(iloc,ival,mode)
    if mode==1
        %===== Filtering with arbitrary resolution =====
        Filter=[ 0  2099 2100 3000 3001 22050 ; 0  0    1    1    0    0 ];
        [syniloc,ind] = sort(iloc);
        FilterEnvelope = interp1(Filter(1,:)',Filter(2,:)',syniloc/N*SR);
        synival = ival(ind)+(20*log10(max(FilterEnvelope,10^-9)));
        synival(ind)=synival;
        syniloc(ind)=syniloc;
    end