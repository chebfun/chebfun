function  [ishappy, epslevel, cutoff] = happinessCheck(f, op, values, vscl, hscl, pref)
%HAPPINESSCHECK   Happiness test for a TRIGTECH
%
% See also CLASSICCHECK, SAMPLETEST.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Grab preferences:
if ( nargin == 1 )
    op = [];
    pref = f.techPref();
elseif ( (nargin == 2) && isstruct(op) )
    pref = op;
    op = [];
elseif ( nargin < 6 )
    pref = f.techPref();
end

if ( nargin < 3 ) 
    values = [];
end
if ( nargin < 4 ) 
    vscl = [];
end
hscl = f.hscale;

% What does happiness mean to you?
if ( strcmpi(pref.happinessCheck, 'classic') )
    % Use the default happiness check procedure from Chebfun V4.
    
    % Check the coefficients are happy:
    [ishappy, epslevel, cutoff] = classicCheck(f, values, vscl, hscl, pref);
    
elseif ( strcmpi(pref.happinessCheck, 'plateau') )
    % Use the 'plateau' happiness check:
    [ishappy, epslevel, cutoff] = plateauCheck(f, values, vscl, hscl, pref);

elseif ( strcmpi(pref.happinessCheck, 'strict') )
    error('CHEBFUN:TRIGTECH:happinessCheck:strictCheck','Strict check not implemented for TRIGTECH.  Please use classic check.');
    
elseif ( strcmpi(pref.happinessCheck, 'loose') )
    error('CHEBFUN:TRIGTECH:happinessCheck:looseCheck','Loose check not implemented for TRIGTECH.  Please use classic check.');
    
else
    % Call a user-defined happiness check:
    [ishappy, epslevel, cutoff] = ...
        pref.happinessCheck(f, values, vscl, hscl,  pref);
    
end

% Check also that sampleTest is happy:
if ( ishappy && ~isempty(op) && ~isnumeric(op) && pref.sampleTest )
    f.epslevel = epslevel;
    ishappy = sampleTest(op, f);
    if ( ~ishappy )
        % It wasn't. Revert cutoff. :(
        cutoff = size(f.values, 1);
    end
end

% set epslevel = eps
epslevel = eps + 0*epslevel;

end
