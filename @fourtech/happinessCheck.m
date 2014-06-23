function  [ishappy, epslevel, cutoff] = happinessCheck(f, op, values, pref)
%HAPPINESSCHECK   Happiness test for a FOURTECH
%
% See also CLASSICCHECK, SAMPLETEST.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Grab preferences:
if ( nargin == 1 )
    op = [];
    pref = f.techPref();
elseif ( (nargin == 2) && isstruct(op) )
    pref = op;
    op = [];
elseif ( nargin < 4 )
    pref = f.techPref();
elseif ( nargin == 3 ) 
    pref = f.techPref(); 
end

% What does happiness mean to you?
if ( strcmpi(pref.happinessCheck, 'classic') )
    % Use the default happiness check procedure from Chebfun V4.
    
    % Check the coefficients are happy:
    [ishappy, epslevel, cutoff] = classicCheck(f, pref);

elseif ( strcmpi(pref.happinessCheck, 'strict') )
    error('CHEBFUN:FOURTECH:happinessCheck:strictCheck','Strict check not implemented for FOURTECH.  Please use classic check.');
    
elseif ( strcmpi(pref.happinessCheck, 'loose') )
    error('CHEBFUN:FOURTECH:happinessCheck:looseCheck','Loose check not implemented for FOURTECH.  Please use classic check.');
    
else
    % Call a user-defined happiness check:
    [ishappy, epslevel, cutoff] = ...
        pref.happinessCheck(f, pref);
    
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

end