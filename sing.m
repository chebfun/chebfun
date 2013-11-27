function song = sing(str,bpm)
% SING - A basic keyboard for MATLAB using chebfuns
%
% S = SING('STR1 STR2 STR3 ...') creates a chebfun F corresponding to the
% musical notes in each of the STR which vaguely correspond to Helmholtz pitch
% notation (http://en.wikipedia.org/wiki/Helmholtz_pitch_notation).
%
% S = SING('STR1 ...', BPM) alters the temp of the tune (in beats per minute).
% The default value is 60bpm.
%
% For example, S = SING('A'), produces a function corresponding to the note A, S
% = SING('B') produces a B, and so on.
%
% Basic chords are provided, for example S = SING('CEG') will produce a major
% triad (http://en.wikipedia.org/wiki/Major_chord). The convention is for the
% kth note to have 1/k times the amplitude of the first.
%
% Empty spaces correspond to pause of a sixteenth note (semiquaver), and commas
% seperate notes without a pause.
%
% If STR is in uppercase, it will last for a quarter note (crotchet), and lower
% case notes for a sixteenth note (semiquaver).
%
% The SING keyboard has three octaves (in A to G). S = SING('C-') produces a low
% C, S = SING('C') produces a middle C, and S = SING('C+') produces a high C.
%
% Sharps are supported by using '#', but there is currently no notation for
% flats.
%
% Examples:
%  1. Ode to Joy:
%     S = SING('BG- BG- CA DB DB CA BG- AD- G-B G-B AD- BG- BG- AD- AD-');
%  2. God Save The Queen
%     S = SING('CG CG DG BG CG DG EC EC FC EC DG CG DG CG BG CG',30);
%
%  See also SOUND, CHEBTUNE.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

s0 = 2^(1/12);
f0 = 440*pi;
fp1 = 2*f0;
fm1 = .5*f0;
if ( nargin < 2 )
    bpm = 60;
    t0 = 1/4;
else
    t0 = 60/bpm/4;
end

q = chebfun('x', [0 t0]);   % quarter note
s = chebfun('x', [0 t0/4]); % sixteenth note

%% middle
% Quarter notes
A = sin(f0*q);
As = sin(f0*s0*q);
B = sin(f0*s0^2*q);
C = sin(f0*s0^3*q);
Cs = sin(f0*s0^4*q);
D = sin(f0*s0^5*q);
Ds = sin(f0*s0^6*q);
E = sin(f0*s0^7*q);
F = sin(f0*s0^8*q);
Fs = sin(f0*s0^9*q);
G = sin(f0*s0^10*q);
Gs = sin(f0*s0^11*q);

% Sixteenth notes
a = sin(f0*s);
as = sin(f0*s0*s);
b = sin(f0*s0^2*s);
c = sin(f0*s0^3*s);
cs = sin(f0*s0^4*s);
d = sin(f0*s0^5*s);
ds = sin(f0*s0^6*s);
e = sin(f0*s0^7*s);
f = sin(f0*s0^8*s);
fs = sin(f0*s0^9*s);
g = sin(f0*s0^10*s);
gs = sin(f0*s0^11*s);

%% -

% Quarter notes
Am = sin(fm1*q);
Ams = sin(fm1*s0*q);
Bm = sin(fm1*s0^2*q);
Cm = sin(fm1*s0^3*q);
Cms = sin(fm1*s0^4*q);
Dm = sin(fm1*s0^5*q);
Dms = sin(fm1*s0^6*q);
Em = sin(fm1*s0^7*q);
Fm = sin(fm1*s0^8*q);
Fms = sin(fm1*s0^9*q);
Gm = sin(fm1*s0^10*q);
Gms = sin(fm1*s0^11*q);

% Sixteenth notes
am = sin(fm1*s);
ams = sin(fm1*s0*s);
bm = sin(fm1*s0^2*s);
cm = sin(fm1*s0^3*s);
cms = sin(fm1*s0^4*s);
dm = sin(fm1*s0^5*s);
dms = sin(fm1*s0^6*s);
em = sin(fm1*s0^7*s);
fm = sin(fm1*s0^8*s);
fm1s = sin(fm1*s0^9*s);
gm = sin(fm1*s0^10*s);
gms = sin(fm1*s0^11*s);

%% +

% Quarter notes
Ap = sin(fp1*q);
Aps = sin(fp1*s0*q);
Bp = sin(fp1*s0^2*q);
Cp = sin(fp1*s0^3*q);
Cps = sin(fp1*s0^4*q);
Dp = sin(fp1*s0^5*q);
Dps = sin(fp1*s0^6*q);
Ep = sin(fp1*s0^7*q);
Fp = sin(fp1*s0^8*q);
Fps = sin(fp1*s0^9*q);
Gp = sin(fp1*s0^10*q);
Gps = sin(fp1*s0^11*q);

% Sixteenth notes
ap = sin(fp1*s);
aps = sin(fp1*s0*s);
bp = sin(fp1*s0^2*s);
cp = sin(fp1*s0^3*s);
cps = sin(fp1*s0^4*s);
dp = sin(fp1*s0^5*s);
dps = sin(fp1*s0^6*s);
ep = sin(fp1*s0^7*s);
fp = sin(fp1*s0^8*s);
fps = sin(fp1*s0^9*s);
gp = sin(fp1*s0^10*s);
gps = sin(fp1*s0^11*s);

%%

if ( nargin == 0 || isnumeric(str) )
    if ( nargin == 0 )
        str = 1;
    end
    switch str
        case 1 % Ode to Joy
            str = 'BG- BG- CA DB DB CA BG- AD- G-B G-B AD- BG- BG- AD- AD-';
        case 2 % God Save The Queen
            str = 'CG CG DG BG CG DG EC EC FC EC DG CG DG CG BG CG';
    end
end

n = length(str);
song = chebfun;
l = 0;
while ( ~isempty(str) )
    l = l+1;
    k = 0;
    for j = 1:numel(str)
        k = k+1;
        if ( strcmp(str(k), ' ') || strcmp(str(k), ',') )
            break
        end
    end
    
    strk = str(1:k);
    str(1:k) = [];
    
    if ( isstrprop(strk(1), 'upper') )
        strk = upper(strk);
        songtmp = 0*A;
    else
        strk = lower(strk);
        songtmp = 0*a;
    end
    
    j = 0;
    while ( ~isempty(strk) )
        if ( length(strk) > 1 && any(strcmpi(strk(2), {'+' '-' '#' '*'})) )
            if ( length(strk) > 2 && any(strcmpi(strk(3), {'#' '*'})) )
                strjk = strk(1:3);
                strk(1:3) = [];
            else
                strjk = strk(1:2);
                strk(1:2) = [];
            end
        else
            strjk = strk(1);
            strk(1) = [];
        end
        j = j+1;
        amp = .4/j;
        
        if ( strcmp(strjk, 'r') || strcmp(strjk, ' ') )
            songtmp = join(songtmp, 0*s);
            break
        end
        
        if strcmp(strjk,'R'), songtmp = join(songtmp, 0*q); break, end
        
        if strcmp(strjk,'A'), songtmp = songtmp + amp*A; end
        if strcmp(strjk,'B'), songtmp = songtmp + amp*B; end
        if strcmp(strjk,'C'), songtmp = songtmp + amp*C; end
        if strcmp(strjk,'D'), songtmp = songtmp + amp*D; end
        if strcmp(strjk,'E'), songtmp = songtmp + amp*E; end
        if strcmp(strjk,'F'), songtmp = songtmp + amp*F; end
        if strcmp(strjk,'G'), songtmp = songtmp + amp*G; end
        
        if strcmp(strjk,'a'), songtmp = songtmp + amp*a; end
        if strcmp(strjk,'b'), songtmp = songtmp + amp*b; end
        if strcmp(strjk,'c'), songtmp = songtmp + amp*c; end
        if strcmp(strjk,'d'), songtmp = songtmp + amp*d; end
        if strcmp(strjk,'e'), songtmp = songtmp + amp*e; end
        if strcmp(strjk,'f'), songtmp = songtmp + amp*f; end
        if strcmp(strjk,'g'), songtmp = songtmp + amp*g; end
        
        if strcmp(strjk,'A-'), songtmp = songtmp + amp*Am; end
        if strcmp(strjk,'B-'), songtmp = songtmp + amp*Bm; end
        if strcmp(strjk,'C-'), songtmp = songtmp + amp*Cm; end
        if strcmp(strjk,'D-'), songtmp = songtmp + amp*Dm; end
        if strcmp(strjk,'E-'), songtmp = songtmp + amp*Em; end
        if strcmp(strjk,'F-'), songtmp = songtmp + amp*Fm; end
        if strcmp(strjk,'G-'), songtmp = songtmp + amp*Gm; end
        
        if strcmp(strjk,'a-'), songtmp = songtmp + amp*am; end
        if strcmp(strjk,'b-'), songtmp = songtmp + amp*bm; end
        if strcmp(strjk,'c-'), songtmp = songtmp + amp*cm; end
        if strcmp(strjk,'d-'), songtmp = songtmp + amp*dm; end
        if strcmp(strjk,'e-'), songtmp = songtmp + amp*em; end
        if strcmp(strjk,'f-'), songtmp = songtmp + amp*fm; end
        if strcmp(strjk,'g-'), songtmp = songtmp + amp*gm; end
        
        if strcmp(strjk,'A+'), songtmp = songtmp + amp*Ap; end
        if strcmp(strjk,'B+'), songtmp = songtmp + amp*Bp; end
        if strcmp(strjk,'C+'), songtmp = songtmp + amp*Cp; end
        if strcmp(strjk,'D+'), songtmp = songtmp + amp*Dp; end
        if strcmp(strjk,'E+'), songtmp = songtmp + amp*Ep; end
        if strcmp(strjk,'F+'), songtmp = songtmp + amp*Fp; end
        if strcmp(strjk,'G+'), songtmp = songtmp + amp*Gp; end
        
        if strcmp(strjk,'a+'), songtmp = songtmp + amp*ap; end
        if strcmp(strjk,'b+'), songtmp = songtmp + amp*bp; end
        if strcmp(strjk,'c+'), songtmp = songtmp + amp*cp; end
        if strcmp(strjk,'d+'), songtmp = songtmp + amp*dp; end
        if strcmp(strjk,'e+'), songtmp = songtmp + amp*ep; end
        if strcmp(strjk,'f+'), songtmp = songtmp + amp*fp; end
        if strcmp(strjk,'g+'), songtmp = songtmp + amp*gp; end
        
        if strcmp(strjk,'a#'), songtmp = songtmp + amp*as; end
        if strcmp(strjk,'c#'), songtmp = songtmp + amp*cs; end
        if strcmp(strjk,'d#'), songtmp = songtmp + amp*ds; end
        if strcmp(strjk,'f#'), songtmp = songtmp + amp*fs; end
        if strcmp(strjk,'g#'), songtmp = songtmp + amp*gs; end
        
        if strcmp(strjk,'A#'), songtmp = songtmp + amp*As; end
        if strcmp(strjk,'C#'), songtmp = songtmp + amp*Cs; end
        if strcmp(strjk,'D#'), songtmp = songtmp + amp*Ds; end
        if strcmp(strjk,'F#'), songtmp = songtmp + amp*Fs; end
        if strcmp(strjk,'G#'), songtmp = songtmp + amp*Gs; end
        
        if strcmp(strjk,'a+#'), songtmp = songtmp + amp*aps; end
        if strcmp(strjk,'c+#'), songtmp = songtmp + amp*cps; end
        if strcmp(strjk,'d+#'), songtmp = songtmp + amp*dps; end
        if strcmp(strjk,'f+#'), songtmp = songtmp + amp*fps; end
        if strcmp(strjk,'g+#'), songtmp = songtmp + amp*gps; end
        
        if strcmp(strjk,'A+#'), songtmp = songtmp + amp*Aps; end
        if strcmp(strjk,'C+#'), songtmp = songtmp + amp*Cps; end
        if strcmp(strjk,'D+#'), songtmp = songtmp + amp*Dps; end
        if strcmp(strjk,'F+#'), songtmp = songtmp + amp*Fps; end
        if strcmp(strjk,'G+#'), songtmp = songtmp + amp*Gps; end
        
        if strcmp(strjk,'a-#'), songtmp = songtmp + amp*ams; end
        if strcmp(strjk,'c-#'), songtmp = songtmp + amp*cms; end
        if strcmp(strjk,'d-#'), songtmp = songtmp + amp*dms; end
        if strcmp(strjk,'f-#'), songtmp = songtmp + amp*fms; end
        if strcmp(strjk,'g-#'), songtmp = songtmp + amp*gms; end
        
        if strcmp(strjk,'A-#'), songtmp = songtmp + amp*Ams; end
        if strcmp(strjk,'C-#'), songtmp = songtmp + amp*Cms; end
        if strcmp(strjk,'D-#'), songtmp = songtmp + amp*Dms; end
        if strcmp(strjk,'F-#'), songtmp = songtmp + amp*Fms; end
        if strcmp(strjk,'G-#'), songtmp = songtmp + amp*Gms; end
        
    end
    song = join(song, songtmp);
end

end
