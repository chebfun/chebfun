function epslevel = updateEpslevel(f, pref)
%UPDATEEPSLEVEL   Determine a new epslevel for a TRIGTECH.
%   UPDATEPSLEVEL(F) calls HAPPINESSCHECK to obtain a new epslevel.
%
%   UPDATEEPSLEVEL(F, EPSBND) assumes that EPSBND is an upper bound for the
%   epslevel of F and passes EPSBND/100 to HAPPINESSCHECK(F, [], [], PREF) as a
%   target epslevel in PREF.EPS.
%
% See also HAPPINESSCHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%   This method is basically called a the end of every operation where we don't
%   know how to update the epslevel manually. This includes, DIFF, CUMSUM, QR,
%   and TIMES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain a preference / epslevel bound:
if ( nargin == 1 )
    pref = f.techPref();
    epslevelBnd = pref.eps;
elseif ( isnumeric(pref) )
    epslevelBnd = pref;
    pref = f.techPref();
    pref.eps = max(epslevelBnd/100, pref.eps);
else
    epslevelBnd = pref.eps;
end

% Call HAPPINESSCHECK()
[ignored, newEpslevel] = happinessCheck(f, [], [], pref);

% Respect the bound:
epslevel = min(newEpslevel, epslevelBnd);

end
    

