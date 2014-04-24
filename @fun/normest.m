function out = normest(f)
%NORMEST   Estimate the Inf-norm of a FUN
%   NORMEST(F) is an estimate of the Inf-norm of the FUN F.

% Call normest of the ONEFUN of F
out = normest(f.onefun);

end
