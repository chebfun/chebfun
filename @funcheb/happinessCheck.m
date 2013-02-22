function  [ishappy, epslevel, cutoff] = happinessCheck(f, op, pref)
%HAPPINESSCHECK Happiness test for a FUNCHEB
%   [ISHAPPY, EPSLEVEL, CUTOFF] = HAPPINESSCHECK(F, OP) tests if the FUNCHEB
%   with values F.VALUES and coefficients F.COEFFS would be a 'happy'
%   approximation (in the sense defined below) to the function handle OP
%   (relative to F.VSCALE and F.HSCALE). If the approximation is happy, the
%   output ISHAPPY is logical TRUE, the happiness level is returned in EPSLEVEL,
%   and CUTOFF indicates the degree to which the coefficients COEFFS may be
%   truncated.
%
%   HAPPINESSCHECK(F, OP, PREF) allows different preferences to be used; in
%   particular PREF.(class(F)).eps as the target tolerance for happiness.
%
%   Furthermore, alternative definitions of happiness can be chosen by setting
%   the PREF.(class(F)).happinessCheck field. This field may be one of the built
%   in checks: 'CLASSIC', 'STRICT', 'LOOSE', or a function handle pointing the a
%   function with the template [ISHAPPY, EPSLEVEL, CUTOFF] = @(F, PREF). The
%   built in check are:
%      CLASSIC: Chooses an EPSLEVEL based upon the length on the 'tail' and a
%               finite difference gradient approximation.
%      STRICT : The tail _must_ be below that specified by PREF.funcheb.eps.
%      LOOSE  : To be specified.
%   Further details of these happiness checks are given in their corresponding
%   m-files.
%
%   Regardless of the happiness definition, HAPPINESSCHECK also performs a
%   SAMPLETEST unless PREF.(class(F).sampleTest is false or OP is empty.
%
% See also classicCheck.m, looseCheck.m, strictCheck.m, sampleTest.m.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Grab preferences:
if ( nargin == 1 )
    op = [];
    pref = f.pref();
elseif ( nargin == 2 && isstruct(op) )
    pref = op;
    op = [];
elseif ( nargin < 3 )
    pref = f.pref();
end

% What does happiness mean to you?
if ( strcmpi(pref.(class(f)).happinessCheck, 'classic') )
    % Use the default happiness check procedure from Chebfun V4.
    
    % Check the coefficients are happy:
    [cutoff, epslevel] = classicCheck(f, pref);
    
    % Happiness here means that the length of the coeffs is decreased:
    if ( cutoff < size(f.values, 1) || cutoff == 1)
        ishappy = true;  % Yay! :)
    else
        ishappy = false; % Boo. :(
    end
    
elseif ( strcmpi(pref.(class(f)).happinessCheck, 'strict') )
    % Use the 'strict' happiness check:
    [ishappy, epslevel, cutoff] = strictCheck(f, pref);
    
elseif ( strcmpi(pref.(class(f)).happinessCheck, 'loose') )
    % Use the 'loose' happiness check:
    [ishappy, epslevel, cutoff] = looseCheck(f, pref);
    
else
    % Call a user-defined happiness check:
    [ishappy, epslevel, cutoff] = ...
        pref.(class(f)).happinessCheck(f, pref);
    
end

% Check also that sampleTest is happy:
if ( ishappy && ~isempty(op) && ~isnumeric(op) && pref.(class(f)).sampletest )
    f.epslevel = epslevel;
    ishappy = funcheb.sampleTest(op, f);
    if ( ~ishappy )
        % It wasn't. Revert cutoff. :(
        cutoff = size(f.values, 1);
    end
end

end