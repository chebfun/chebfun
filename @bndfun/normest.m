function nm = normest(f)
% NORMEST Estimate the Inf-norm of a BNDFUN
%    NORMEST(F) is an estimate of the Inf-norm of the BNDFUN F.

% Call normest of the ONEFUN of F
nm = normest(f.onefun);

end