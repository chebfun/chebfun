function  [ishappy, epslevel, cutoff] = happinessCheck(f, op, pref)
%HAPPINESSCHECK   Happiness test for a CHEBTECH
%   [ISHAPPY, EPSLEVEL, CUTOFF] = HAPPINESSCHECK(F, OP) tests if the CHEBTECH
%   with values F.VALUES and coefficients F.COEFFS would be a 'happy'
%   approximation (in the sense defined below and relative to F.VSCALE and
%   F.HSCALE) to the function handle OP. If the approximation is happy, the
%   output ISHAPPY is TRUE, the happiness level is returned in EPSLEVEL,
%   and CUTOFF indicates the point to which the coefficients COEFFS may be
%   truncated. Even if ISHAPPY is FALSE, the attempted happiness level is still
%   returned in EPSLEVEL (i.e., we attempted to be happy at EPSLEVEL but failed)
%   and CUTOFF is returned as size(f.values, 1).
%
%   HAPPINESSCHECK(F, OP, PREF) allows different preferences to be used; in
%   particular PREF.CHEBTECH.EPS sets the target tolerance for happiness.
%
%   Furthermore, alternative definitions of happiness can be chosen by setting
%   the PREF.CHEBTECH.HAPPINESSCHECK field. This field may be one of the built
%   in checks: 'CLASSIC', 'STRICT', 'LOOSE', or a function handle pointing to a
%   function with the template [ISHAPPY, EPSLEVEL, CUTOFF] = @(F, PREF). The
%   built in checks are:
%      CLASSIC: Chooses an EPSLEVEL based upon the length on the 'tail' and a
%               finite difference gradient approximation.
%      STRICT : The tail _must_ be below that specified by PREF.CHEBTECH.EPS.
%      LOOSE  : To be specified.
%   Further details of these happiness checks are given in their corresponding
%   m-files.
%
%   Regardless of the happiness definition, HAPPINESSCHECK also performs a
%   SAMPLETEST unless PREF.CHEBTECH.SAMPLETEST is FALSE or OP is empty.
%
% See also CLASSICCHECK, LOOSECHECK, STRICTCHECK, SAMPLETEST.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Grab preferences:
if ( nargin == 1 )
    op = [];
    pref = f.pref();
elseif ( (nargin == 2) && isstruct(op) )
    pref = op;
    op = [];
elseif ( nargin < 3 )
    pref = f.pref();
end

% What does happiness mean to you?
if ( strcmpi(pref.chebtech.happinessCheck, 'classic') )
    % Use the default happiness check procedure from Chebfun V4.
    
    % Check the coefficients are happy:
    [ishappy, epslevel, cutoff] = classicCheck(f, pref);

elseif ( strcmpi(pref.chebtech.happinessCheck, 'strict') )
    % Use the 'strict' happiness check:
    [ishappy, epslevel, cutoff] = strictCheck(f, pref);
    
elseif ( strcmpi(pref.chebtech.happinessCheck, 'loose') )
    % Use the 'loose' happiness check:
    [ishappy, epslevel, cutoff] = looseCheck(f, pref);
    
else
    % Call a user-defined happiness check:
    [ishappy, epslevel, cutoff] = ...
        pref.chebtech.happinessCheck(f, pref);
    
end

% Check also that sampleTest is happy:
if ( ishappy && ~isempty(op) && ~isnumeric(op) && pref.chebtech.sampletest )
    f.epslevel = epslevel;
    ishappy = sampleTest(op, f);
    if ( ~ishappy )
        % It wasn't. Revert cutoff. :(
        cutoff = size(f.values, 1);
    end
end

end
