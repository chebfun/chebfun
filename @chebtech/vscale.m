function vscl = vscale(f)
%VSCALE  Vertical scale of a CHEBTECH. 
%   VSCALE(F) returns a row vector storing the magnitude of the largest entry in
%   each column of the CHEBTECH sampled on its Chebyshev grid.

values = f.coeffs2vals(f.coeffs);
vscl = max(abs(values), [], 1);

end