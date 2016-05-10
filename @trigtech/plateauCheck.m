function [ishappy, cutoff] = plateauCheck(f, values, data, pref)
%PLATEAUCHECK   Attempt to trim trailing TRIGIER coefficients in a TRIGTECH.
%   [ISHAPPY, CUTOFF] = PLATEAUCHECK(F, VALUES, DATA) returns an estimated
%   location, the CUTOFF, at which the TRIGTECH F could be truncated. One of
%   two criteria must be met: Either:
%
%     (1) The coefficients are sufficiently small (as specified by the default
%     EPS property of TRIGTECH) relative to F.VSCALE (or using absolute size if
%     VSCALE(F)=0); or
%
%     (2) The coefficients are somewhat small and apparently unlikely to
%     continue decreasing in a meaningful amount (i.e., have reached a "plateau"
%     in convergence).
%
%   The reason for criterion (2) is that the problem may have a large condition
%   number that prevents convergence to the full requested accuracy, as often
%   happens in collocation of differential equations.
%
%   Output CUTOFF is an estimate of how many of the coefficients are useful.
%
%   [ISHAPPY, CUTOFF] = PLATEAUCHECK(F, VALUES, DATA, PREF) allows additional
%   preferences to be passed. In particular, one can adjust the target accuracy
%   with PREF.CHEBFUNEPS.
%
% See also LINOPV4CHECK, STRICTCHECK, CLASSICCHECK.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: implement PLATEAUCHECK for TRIGTECH. For the moment, we just call
% classicCheck. The reason why the plateauCheck() is needed for TRIGTECH is 
% that it gets called in OPDISCRETIZATION/TESTCONVERGENCE.

[ishappy, cutoff] = classicCheck(f, values, data, pref);

end
