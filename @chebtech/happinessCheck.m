function  [ishappy, epslevel, cutoff] = happinessCheck(f, op, values, vscl, hscl, pref)
%HAPPINESSCHECK   Happiness test for a CHEBTECH
%   [ISHAPPY, EPSLEVEL, CUTOFF] = HAPPINESSCHECK(F, OP, VALUES, VSCL, HSCL)
%   tests if the CHEBTECH with values VALUES and coefficients F.COEFFS would be
%   a 'happy' approximation (in the sense defined below and relative to VSCL and
%   HSCL) to the function handle OP. If the approximation is happy, the output
%   ISHAPPY is TRUE, the happiness level is returned in EPSLEVEL, and CUTOFF
%   indicates the point to which the coefficients COEFFS may be truncated. If
%   ISHAPPY is false, EPSLEVEL returns an estimate of the accuracy achieved.
%
%   HAPPINESSCHECK(F) computes VALUES used above from F.COEFFS2VALS(F.COEFFS).
%
%   HAPPINESSCHECK(..., PREF) allows different preferences to be used; in
%   particular PREF.EPS sets the target tolerance for happiness.  If
%   constructing an array-valued CHEBTECH, PREF.EPS may be a row vector of
%   target tolerances for each column.
%
%   Furthermore, alternative definitions of happiness can be chosen by setting
%   the PREF.HAPPINESSCHECK field. This field may be one of the built in
%   checks: 'CLASSIC', 'STRICT', 'LOOSE', or a function handle pointing to a
%   function with the template [ISHAPPY, EPSLEVEL, CUTOFF] = @(F, PREF). The
%   built in checks are:
%      CLASSIC: Chooses an EPSLEVEL based upon the length on the 'tail' and a
%               finite difference gradient approximation.
%      STRICT : The tail _must_ be below that specified by PREF.EPS.
%      LOOSE  : To be specified.
%   Further details of these happiness checks are given in their corresponding
%   m-files.
%
%   Regardless of the happiness definition, HAPPINESSCHECK also performs a
%   SAMPLETEST unless PREF.SAMPLETEST is FALSE or OP is empty.
%
% See also CLASSICCHECK, LOOSECHECK, STRICTCHECK, SAMPLETEST.

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
%     values = f.coeffs2vals(f.coeffs);
end
if ( nargin < 4 )
    vscl = [];
end
if ( nargin < 4 )
    hscl = [];
end

% What does happiness mean to you?
if ( strcmpi(pref.happinessCheck, 'classic') )
    % Use the default happiness check procedure from Chebfun V4.
    
    % Check the coefficients are happy:
    [ishappy, epslevel, cutoff] = classicCheck(f, values, vscl, hscl, pref);

elseif ( strcmpi(pref.happinessCheck, 'strict') )
    % Use the 'strict' happiness check:
    [ishappy, epslevel, cutoff] = strictCheck(f, values, vscl, hscl, pref);
    
elseif ( strcmpi(pref.happinessCheck, 'loose') )
    % Use the 'loose' happiness check:
    [ishappy, epslevel, cutoff] = looseCheck(f, values, vscl, hscl, pref);
    
elseif ( strcmpi(pref.happinessCheck, 'plateau') )
    % Use the 'plateau' happiness check:
    [ishappy, epslevel, cutoff] = plateauCheck(f, values, vscl, hscl, pref);    
    
else
    % Call a user-defined happiness check:
    [ishappy, epslevel, cutoff] = ...
        pref.happinessCheck(f, values, vscl, hscl, pref);
    
end

% Check also that sampleTest is happy:
if ( ishappy && ~isempty(op) && ~isnumeric(op) && pref.sampleTest )
    f.epslevel = epslevel;
    ishappy = sampleTest(op, values, f, vscl);
    if ( ~ishappy )
        % It wasn't. Revert cutoff. :(
        cutoff = size(values, 1);
    end
end

% set epslevel = eps
epslevel = eps + 0*epslevel;

end
