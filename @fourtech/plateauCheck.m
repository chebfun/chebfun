function [ishappy, epsLevel, cutoff] = plateauCheck(f, pref)
%PLATEAUCHECK   Attempt to trim trailing FOURIER coefficients in a FOURTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = PLATEAUCHECK(F, VALUES) returns an estimated
%   location, the CUTOFF, at which the FOURTECH F could be truncated. One of two
%   criteria must be met: Either:
%
%     (1) The coefficients are sufficiently small (as specified by the default
%     EPS property of FOURTECH) relative to F.VSCALE (or using absolute size if
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: implement PLATEAUCHECK for FOURTECH.
[ishappy, epsLevel, cutoff] = classicCheck(f, pref);

end