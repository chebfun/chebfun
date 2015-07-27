function [ishappy, epslevel, cutoff] = plateauCheck(f, values, vscl, hscl, pref)
%PLATEAUCHECK   Attempt to trim trailing TRIGIER coefficients in a TRIGTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = PLATEAUCHECK(F, VALUES) returns an estimated
%   location, the CUTOFF, at which the TRIGTECH F could be truncated. One of two
%   criteria must be met: Either:
%
%     (1) The coefficients are sufficiently small (as specified by the default
%     EPS property of TRIGTECH) relative to F.VSCALE (or using absolute size if
%     F.VSCALE=0); or
%
%     (2) The coefficients are somewhat small and apparently unlikely to
%     continue decreasing in a meaningful amount (i.e., have reached a "plateau"
%     in convergence).
%
%   The reason for criterion (2) is that the problem may have a large condition
%   number that prevents convergence to the full requested accuracy, as often
%   happens in the collocation of differential equations.
%
%   Output EPSLEVEL is an estimate of the relative size of the last
%   "meaningful" expansion coefficients of the function, and the output 
%   CUTOFF is an estimate of how many of the coefficients are useful.
%
%   [ISHAPPY, EPSLEVEL, CUTOFF] = PLATEAUCHECK(F, VALUES, PREF) allows
%   additional preferences to be passed. In particular, one can adjust the
%   target accuracy with PREF.EPS.
%
% See also LINOPV4CHECK, STRICTCHECK, CLASSICCHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: We currently need to do this as trigtech has a vscale and chebtech
% doesn't.
f.vscale = max(f.vscale, vscl);

% [TODO]: implement PLATEAUCHECK for TRIGTECH. For the moment, we just call
% classicCheck. The reason why the plateauCheck() is needed for TRIGTECH is 
% that it gets called in CHEBDISCRETIZATION/TESTCONVERGENCE.

[ishappy, epslevel, cutoff] = classicCheck(f, values, vscl, f.hscale, pref);

end
