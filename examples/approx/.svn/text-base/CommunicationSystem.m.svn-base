%% Illustrating the mathematics of signal processing in Chebfun
% Mohsin Javed, 24th August, 2012

%%
% (Chebfun example approx/CommunicationSystem.m)
% [Tags: #signaprocessing, #AM, #amplitudemodulation, #DIRAC, #deltafunction]

%% Introduction
% In this example we use Chebfun to illustrate some of the basic
% mathematics of communication systems, including the interpretation of
% discrete signals as combinations of Dirac delta functions.  We hasten to
% add that Chebfun is absolutely not a practical tool for signal
% processing! It provides a convenient setting, however, for the
% exploration of basic principles.
% 
% Typically the purpose of a communication system is to transmit a low pass
% signal form one location to another through a medium which allows
% efficient propagation of electromagnetic waves. The medium is usually
% referred to as the `channel' and it may be a telephone land-line, a fiber
% optic cable or air etc. It is well known that to transmit a signal
% efficiently, it must be modulated at the transmitter end with a carrier
% wave of high frequency [1]. At the receiver end the modulated signal is
% then `demodulated' to recover the original signal. The demodulation
% process typically involves low pass filtering to reconstruct the original
% signal.

%% Amplitude Modulation
% There are various types of modulation/demodulation schemes but the
% simplest one is known as amplitude modulation or AM. In AM, the original
% information is carried in the amplitude of the carrier wave. Let us try
% to understand the basics of AM using Chebfun.
%
% Throughout this example, we will represent functions in the time domain
% with small letters and functions in the frequency domain with capital
% letters.

%% A Low Frequency Signal
% We start with the construction of a low frequency signal with maximum
% angular frequency, say
wmax = 4;

%%
% We now select some random Fourier coefficients

a  = -1+2*rand(wmax,1);   % Amplitude of cosines
b  = -1+2*rand(wmax,1);   % Amplitude of sines
a0 = rand;                % The DC-component

%%
% We can use these coefficients to construct
% a low frequency random signal as a chebfun.
d = domain(-5,5);
t = chebfun('t',d);
s = chebfun(a0,d);
for k = 1:wmax
    s = s + a(k)*cos(k*t)+b(k)*sin(k*t);
end
plot(s)
%%
% Perhaps it is interesting to note that we can not only see but actually
% hear this chebfun as well.
chebtune(s)

%% Fourier Space
% Before we make any attempts to transmit this signal, let us view the
% signal in the frequency domain. This means that we need to find the
% Fourier transform of the signal. Unfortunately, Chebfun can not find the
% Fourier transform of a signal but for this tailor-made signal, it is not
% to hard to construct its Fourier transform. First of all, we recall that
% the Fourier transform of a constant function is the Dirac-delta function.
wd = domain(-5*wmax-1,5*wmax+1);
w = chebfun('w',wd);
S = a0*dirac(w);

%%
% The Fourier transform of cosines (real and even) and sines (real
% and odd) come in pairs of real/even and imaginary/odd
% Dirac-impulses respectively.
Fcos = @(k) (dirac(w-k) + dirac(w+k))/2;
Fsin = @(k) (dirac(w-k) - dirac(w+k))/2i;

%%
% We can now construct the Fourier transform of the signal easily.
for k = 1:wmax
    S = S + a(k)*Fcos(k)+b(k)*Fsin(k);
end

%%
% Let us see how it looks. We first plot the frequency components
% corresponding to the even part of the signal.
plot(real(S)), xlim([wd(1), wd(2)])
title('Spectrum of the even part of the signal')

%%
% And then plot the frequency components corresponding to the odd part of
% the signal.
plot(imag(S),'r'); xlim([wd(1), wd(2)])
title('Spectrum of the odd part of the signal')

%%
% It is easy to see that all the blue arrows (cosines) are in pairs,
% pointing in the same direction, while pairs of red arrows (sines) are
% pointing in the opposite direction. The blue arrow in the centre is the
% DC component of the signal.

%% Modulation 
% We now modulate this signal with a high frequency carrier wave. In AM,
% this is accomplished by simply multiplying the signal with a high
% frequency sinusoidal signal.
wc = 5*wmax;             % carrier frequency
stx = s.*cos(wc*t);      % transmitted signal
plot(stx); hold on

%%
% This multiplication has the effect of 'chopping' the original signal
% several times. We can actually see the original signal as the envelope of
% this modulated signal, hence the name amplitude modulation.
plot(s, '-.r'), hold off

%% 
% And one might be interested in listening to the modulated signal as well.
chebtune(stx)

%%
% Let's see how the spectrum of the modulated signal has changed due to
% this multiplication. The convolution theorem tells us that multiplication
% in time domain is convolution in frequency domain.
Stw = conv(S,Fcos(wc));
plot(real(Stw)), xlim([Stw.ends(1) Stw.ends(end)]), hold on
plot(imag(Stw),'r'), xlim([Stw.ends(1) Stw.ends(end)]), hold off
title('Spectrum of the modulated signal (real and imaginary part on the same axis)')

%%
% We can see that the spectrum is now nicely centered around the carrier
% frequency and is now ready for transmission. Once the transmission is
% done, the electromagnetic waves happily propagate with lightening speed,
% quite literally, towards their destination.

%% Demodulation
% At the receiver end, we now try to demodulate the signal. In AM, this is
% simply done my multiplying the received signal again with a sinusoidal of
% the same frequency as the one used for modulation.
srx = stx.*cos(wc*t);
plot(srx); hold on
plot(s,'-.r'), hold off

%%
% As usual, let's see what is happening in the frequency domain. The
% demodulation corresponds to another convolution with the Fourier
% transform of the carrier wave.
Srw = conv(Stw,Fcos(wc));
plot(real(Srw)), xlim([Srw.ends(1) Srw.ends(end)]), hold on
plot(imag(Srw),'r'), xlim([Srw.ends(1) Srw.ends(end)]), hold off
title('Spectrum of the demodulated signal (real and imaginary part on the same axis)')


%% Low Pass Filtering
% To recover the original signal, we now need to get rid of the high
% frequency components shown in the above two plots. There are various ways
% to do this but all of them involve some kind of low pass filtering.
% Mathematically, we can do the same by convolving the received signal with
% the impulse response of a low pass filter. Let us use the the ideal low
% pass filter, whose impulse response is well known to be the sinc
% function. To make sure that we get all the frequencies of the original
% signal, we use a slightly greater cutoff frequency in the argument of the
% sinc filter:
sinc = sin((wmax+2)*t)./(pi*t);

%%
% We now apply this filter to the received signal.
sf = 2*conv(srx,sinc);

%%
% One might wonder why a factor of 2 has been used in the convolution
% above. The answer becomes apparent if one recalls that the original
% signal has effectively been multiplied with the square of a cosine and $$
% \cos^2 t = \frac{(1+\cos 2t )}{2}. $$
% 
%% 
% Since convolution in Chebfun extends the domain of the signal, we only
% retain the interval corresponding to the original domain.
sf = sf{d(1),d(2)};
plot(sf), hold on

%%
% How well does this compare with the original signal in the eyeball-norm?
plot(s,'-.r'), hold off

%%
% The match seems good in the middle and poor at the ends. This was
% expected because the ideal filtering process is done on an infinite
% domain with no aliasing.

%%
% How does this compare with the original signal in the ear-norm? Let's
% listen to both of them. To refresh your memory, here is the original
% signal once again.
chebtune(s)

%%
% And here goes the received signal.
chebtune(sf)


%%
% Finally, we illustrate the filtering process in the frequency domain as
% well. Recall that the Fourier transform of the sinc function is a
% rectangular pulse. Therefore low pass filtering corresponds to a
% multiplication by this pulse in the frequency domain.
rectw = @(w) heaviside(w+wmax+2)-heaviside(w-wmax-2);
rect = chebfun(@(w) rectw(w), domain(Srw),'splitting', 'on');
Sf = rect.*Srw;
plot(real(Sf)), hold on
plot(imag(Sf), 'r')
title('Spectrum of the filtered signal (real and imaginary part on the same axis)')
plot(rect,'k', 'jumpline', '-')
ylim([-1.2 1.2]), xlim([Sf.ends(1) Sf.ends(end)]), hold off

%%
% We again remind the reader that this example does not suggest using
% Chebfun for signal processing. However, the aim here is to illustrate the
% ease with which functional computing can make us understand mathematical
% concepts in a variety of practical problems.

%%
% References
%
% [1] Oppenheim and Willsky, Signals and Systems, second edition, Prentice
% Hall, 1998.