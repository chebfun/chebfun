function song = sing(str, bpm)
%SING   A basic keyboard for MATLAB using CHEBFUNs.
%   S = SING('STR1 STR2 STR3 ...') creates a CHEBFUN F corresponding to the
%   musical notes in each of the STR which vaguely correspond to Helmholtz pitch
%   notation (http://en.wikipedia.org/wiki/Helmholtz_pitch_notation).
%
%   S = SING('STR1 ...', BPM) alters the temp of the tune (in beats per minute).
%   The default value is 60bpm.
%
%   For example, S = SING('A'), produces a function corresponding to the note
%   A, S = SING('B') produces a B, and so on.
%
%   Basic chords are provided, for example S = SING('CEG') will produce a major
%   triad (http://en.wikipedia.org/wiki/Major_chord). The convention is for the
%   kth note to have 1/k times the amplitude of the first.
%
%   Empty spaces correspond to pause of a sixteenth note (semiquaver), and
%   commas separate notes without a pause.
%
%   If STR is in uppercase, it will last for a quarter note (crotchet), and
%   lower case notes for a sixteenth note (semiquaver).
%
%   The SING keyboard has three octaves (in A to G). S = SING('C-') produces a
%   low C, S = SING('C') produces a middle C, and S = SING('C+') produces a
%   high C.
%
%   Sharps are supported by using '#', but there is currently no notation for
%   flats.
%
%   Examples:
%    1. Ode to Joy:
%       S = SING('BG- BG- CA DB DB CA BG- AD- G-B G-B AD- BG- BG- AD- AD-');
%    2. God Save The Queen
%       S = SING('CG CG DG BG CG DG EC EC FC EC DG CG DG CG BG CG',30);
%
%   CHEBFUNs for these examples may be generated with SING(1) and SING(2),
%   respectively.
%
% See also SOUND, CHEBTUNE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set up basic frequencies.
s0 = 2^(1/12);
f0 = 440*pi;
fp1 = 2*f0;
fm1 = .5*f0;

% Set tempo.
if ( nargin < 2 )
    bpm = 60;
    t0 = 1/4;
else
    t0 = 60/bpm/4;
end

q = chebfun('x', [0 t0]);   % Quarter note.
s = chebfun('x', [0 t0/4]); % Sixteenth note.

%% Middle

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
fms = sin(fm1*s0^9*s);
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

% Easy access to example tunes.
if ( (nargin == 0) || isnumeric(str) )
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

% Build the CHEBFUN from the string of notes for the tune.
song = chebfun;
l = 0;
while ( ~isempty(str) )
    l = l + 1;
    k = 0;

    % Get the string for the next chord.
    for j = 1:numel(str)
        k = k + 1;
        if ( strcmp(str(k), ' ') || strcmp(str(k), ',') )
            break
        end
    end

    strk = str(1:k);
    str(1:k) = [];

    % Quarter note or sixteenth note?
    if ( isstrprop(strk(1), 'upper') )
        strk = upper(strk);
        songtmp = 0*A;
    else
        strk = lower(strk);
        songtmp = 0*a;
    end

    % Assemble the chord from its constituent notes.
    j = 0;
    while ( ~isempty(strk) )
        % Parse out one note from the chord.
        if ( (length(strk) > 1) && any(strcmpi(strk(2), {'+' '-' '#' '*'})) )
            if ( (length(strk) > 2) && any(strcmpi(strk(3), {'#' '*'})) )
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

        % Set the amplitude of this note in the chord.
        j = j + 1;
        amp = .4/j;

        % Decode and add in the note.
        if ( strcmp(strjk, 'r') || strcmp(strjk, ' ') )  % Sixteenth rest.
            songtmp = join(songtmp, 0*s);
            break
        elseif ( strcmp(strjk, 'R') )                    % Quarter rest.
            songtmp = join(songtmp, 0*q);
            break
        else                                             % Actual note.
            songtmp = songtmp + amp*getNote(strjk);
        end
    end

    % Tack the chord onto the end of the song.
    song = join(song, songtmp);
end

    % Nested function for mapping note strings to note CHEBFUNs defined above.
    function note = getNote(s)
        switch (s)
            case 'A', note = A;
            case 'B', note = B;
            case 'C', note = C;
            case 'D', note = D;
            case 'E', note = E;
            case 'F', note = F;
            case 'G', note = G;

            case 'a', note = a;
            case 'b', note = b;
            case 'c', note = c;
            case 'd', note = d;
            case 'e', note = e;
            case 'f', note = f;
            case 'g', note = g;

            case 'A-', note = Am;
            case 'B-', note = Bm;
            case 'C-', note = Cm;
            case 'D-', note = Dm;
            case 'E-', note = Em;
            case 'F-', note = Fm;
            case 'G-', note = Gm;

            case 'a-', note = am;
            case 'b-', note = bm;
            case 'c-', note = cm;
            case 'd-', note = dm;
            case 'e-', note = em;
            case 'f-', note = fm;
            case 'g-', note = gm;

            case 'A+', note = Ap;
            case 'B+', note = Bp;
            case 'C+', note = Cp;
            case 'D+', note = Dp;
            case 'E+', note = Ep;
            case 'F+', note = Fp;
            case 'G+', note = Gp;

            case 'a+', note = ap;
            case 'b+', note = bp;
            case 'c+', note = cp;
            case 'd+', note = dp;
            case 'e+', note = ep;
            case 'f+', note = fp;
            case 'g+', note = gp;

            case 'a#', note = as;
            case 'c#', note = cs;
            case 'd#', note = ds;
            case 'f#', note = fs;
            case 'g#', note = gs;

            case 'A#', note = As;
            case 'C#', note = Cs;
            case 'D#', note = Ds;
            case 'F#', note = Fs;
            case 'G#', note = Gs;

            case 'a+#', note = aps;
            case 'c+#', note = cps;
            case 'd+#', note = dps;
            case 'f+#', note = fps;
            case 'g+#', note = gps;

            case 'A+#', note = Aps;
            case 'C+#', note = Cps;
            case 'D+#', note = Dps;
            case 'F+#', note = Fps;
            case 'G+#', note = Gps;

            case 'a-#', note = ams;
            case 'c-#', note = cms;
            case 'd-#', note = dms;
            case 'f-#', note = fms;
            case 'g-#', note = gms;

            case 'A-#', note = Ams;
            case 'C-#', note = Cms;
            case 'D-#', note = Dms;
            case 'F-#', note = Fms;
            case 'G-#', note = Gms;
        end
    end

end
