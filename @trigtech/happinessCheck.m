function  [ishappy, cutoff] = happinessCheck(f, op, values, data, pref)
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
elseif ( nargin < 5 )
    pref = f.techPref();
end

if ( nargin < 3 ) 
    values = [];
end
if ( nargin < 4 ) 
    data = struct();
end
data = trigtech.parseDataInputs(data);

% vscale defaults to zero if not given but should be at least f.vscale.  (Only
% makes sense to have a larger "global" vscale.)
data.vscale = max(data.vscale, vscale(f));

% What does happiness mean to you?
if ( strcmpi(pref.happinessCheck, 'standard') )
    % Use the 'standard' happiness check:
    [ishappy, cutoff] = standardCheck(f, values, data, pref);

elseif ( strcmpi(pref.happinessCheck, 'classic') )
    % Use the default happiness check procedure from Chebfun V4.
    [ishappy, cutoff] = classicCheck(f, values, data, pref);
    
elseif ( strcmpi(pref.happinessCheck, 'plateau') )
    % Use the 'plateau' happiness check:
    [ishappy, cutoff] = plateauCheck(f, values, data, pref);

elseif ( strcmpi(pref.happinessCheck, 'strict') )
    error('CHEBFUN:TRIGTECH:happinessCheck:strictCheck',...
          'Strict check not implemented for TRIGTECH.  Please use classic check.');
    
elseif ( strcmpi(pref.happinessCheck, 'loose') )
    error('CHEBFUN:TRIGTECH:happinessCheck:looseCheck',...
           'Loose check not implemented for TRIGTECH.  Please use classic check.');
    
else
    % Call a user-defined happiness check:
    checker = pref.happinessCheck;
    [ishappy, cutoff] = checker(f, values, data, pref);
    
end

% Check also that sampleTest is happy:
if ( ishappy && ~isempty(op) && ~isnumeric(op) && pref.sampleTest )
    ishappy = sampleTest(op, f, pref);
    if ( ~ishappy )
        % It wasn't. Revert cutoff. :(
        cutoff = size(f.values, 1);
    end
end

end
