function varargout = chebtune(f, d)
%CHEBTUNE   CHEBFUN melody player.
%   CHEBTUNE(F) plays melodies with varying pitch corresponding to the real part
%   of the function values of each CHEBFUN in F. The function value 0 is
%   associated with the tone c'' and the integers below and above correspond to
%   the semi-tones. The melodies are separated in the stereo panorama.
%
%   CHEBTUNE(F, D) plays the melody for D seconds. The default value is D = 2.
%
% Example: CHEBPOLY-phony
%      f = 7*chebpoly(0:2) - 7;
%      f = [f , f + .2];  % add some chorus
%      chebtune(f, 3);
%
% Example: Police car
%      x = chebfun('x');
%      chebtune([9 + 6*sin(46*x), 7 + 10*sin(4*x)], 5);
%
% Example: Can you hear the shape of a CHEBFUN?
%      f = chebfun(12*rand(6, 1) - 6);
%      chebtune(f, 6);
%      plot(f);
%
% See also SOUND, SING.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    d = 2;
end

d = max(d, .5);
Fc = 22050;          % Sound sampling rate.
df = 75;             % CHEBFUN sampling divisor.
dom = f.domain;

t = linspace(dom(1), dom(end), d*Fc/df).';
s = real(feval(f, t));             % Sample
s(abs(s) > 60 | isnan(s)) = -Inf;  % Silence
tone0 = 523.25;                    % c''
s = tone0*2.^(s/12);               % Note to frequency.
s = kron(s, ones(df, 1));          % Upsampling.
s = .5*sin(2*pi*cumsum(s)/Fc);     % Frequency to waveform.
ch = size(s, 2);                   % Number of channels.

% Fade in and out to remove blub:
fl = 500;
fade = repmat(linspace(1, 0, fl).', 1, ch);
s(1:fl,:) = s(1:fl,:).*fade(end:-1:1,:);
s(end+1-fl:end,:) = s(end+1-fl:end,:).*fade;
s = [s ; zeros(5000, ch)];
if ( ch > 1 )
    chv = linspace(0.1, .9, ch).'; 
    chv = chv/sum(chv);
    s = s*[chv, chv(end:-1:1)];
else
    s = [s , s];
end

% Use Java audio player to play back our sound.  (sound(s, Fc) playblocks.)
if ( usejava('jvm') )
    persistent ap
    ap = audioplayer(s, Fc);
    play(ap);
else
    warning('CHEBFUN:CHEBFUN:chebtune:unsupportedPlatform', ...
        ['This platform does not support specifing FS or BITS when not ', ...
         'using Java.']);
end

if ( nargout > 0)
    varargout(1) = {s}; 
end

if ( nargout > 1 ) 
    varargout(2) = {Fc}; 
end

end
