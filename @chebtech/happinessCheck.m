function  [ishappy, cutoff] = happinessCheck(f, op, values, data, pref)
%HAPPINESSCHECK   Happiness test for a CHEBTECH
%   [ISHAPPY, CUTOFF] = HAPPINESSCHECK(F, OP, VALUES, DATA) tests if the
%   CHEBTECH with values VALUES and coefficients F.COEFFS would be a 'happy'
%   approximation (in the sense defined below and relative to DATA.VSCALE 
%   and DATA.HSCALE) to the function handle OP. If the approximation is 
%   happy, the output ISHAPPY is TRUE, the  CUTOFF indicates the point to 
%   which the coefficients COEFFS may be truncated.
%
%   HAPPINESSCHECK(F) computes VALUES used above from F.COEFFS2VALS(F.COEFFS).
%
%   HAPPINESSCHECK(F, OP, VALUES, DATA,PREF) allows different preferences 
%   to be used; in particular PREF.CHEBFUNEPS sets the target tolerance for
%   happiness.  If constructing an array-valued CHEBTECH, PREF.CHEBFUNEPS 
%   may be a row vector of target tolerances for each column.
%
%   Furthermore, alternative definitions of happiness can be chosen by setting
%   the PREF.HAPPINESSCHECK field. This field may be one of the built in
%   checks: 'CLASSIC', 'STRICT', 'LOOSE', or a function handle pointing to a
%   function with the template [ISHAPPY, CUTOFF] = @(F, PREF). The built in
%   checks are:
%      CLASSIC: Chooses an CUTOFF based upon the length on the 'tail' and a
%               finite difference gradient approximation.
%      STRICT : The tail _must_ be below that specified by PREF.CHEBFUNEPS.
%      LOOSE  : To be specified.
%   Further details of these happiness checks are given in their corresponding
%   m-files.
%
%   Regardless of the happiness definition, HAPPINESSCHECK also performs a
%   SAMPLETEST unless PREF.SAMPLETEST is FALSE or OP is empty.
%
% See also CLASSICCHECK, LOOSECHECK, STRICTCHECK, SAMPLETEST.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
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
data = chebtech.parseDataInputs(data, pref);

% vscale defaults to zero if not given but should be at least vscale(f). 
% (It only makes sense to have a larger "global" vscale.)
data.vscale = max(data.vscale, vscale(f));

% What does happiness mean to you?
if ( strcmpi(pref.happinessCheck, 'standard') )
    % Use the 'standard' happiness check:
    [ishappy, cutoff] = standardCheck(f, values, data, pref);

elseif ( strcmpi(pref.happinessCheck, 'classic') )
    % Use the default happiness check procedure from Chebfun V4.
    [ishappy, cutoff] = classicCheck(f, values, data, pref);

elseif ( strcmpi(pref.happinessCheck, 'strict') )
    % Use the 'strict' happiness check:
    [ishappy, cutoff] = strictCheck(f, values, data, pref);
    
elseif ( strcmpi(pref.happinessCheck, 'loose') )
    % Use the 'loose' happiness check:
    [ishappy, cutoff] = looseCheck(f, values, data, pref);
    
elseif ( strcmpi(pref.happinessCheck, 'plateau') )
    % Use the 'plateau' happiness check:
    [ishappy, cutoff] = plateauCheck(f, values, data, pref);
    
else
    % Call a user-defined happiness check:
    [ishappy, cutoff] = pref.happinessCheck(f, values, data, pref);
    
end

% Check also that sampleTest is happy:
if ( ishappy && ~isempty(op) && ~isnumeric(op) && pref.sampleTest )
    ishappy = sampleTest(op, values, f, data, pref);
    if ( ~ishappy )
        % It wasn't. Revert cutoff. :(
        cutoff = size(values, 1);
    end
end

end
