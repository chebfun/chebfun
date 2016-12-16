function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Creates plotcoeffs of a CHEBFUN3T object.
%   PLOTCOEFFS(F) creates a scatter3 plot of the tensor of coefficients of 
%   the CHEBFUN3T object F that visualizes logarithm of the magnitude of 
%   each entry.
%
% See also CHEBFUN3/PLOTCOEFFS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% The coefficients and vertical scale:
absCoeffs = abs(f.coeffs);
vscl = f.vscale;

% Get the size:
[m,n,p] = size(absCoeffs);

% Generate vectors needed for the scatter3 plot
[ind1, ind2, ind3] = ndgrid(1:m, 1:n, 1:p);
x = ind1(:);
y = ind2(:); 
z = ind3(:);

% Add a tiny amount to zeros to make plots look nicer:
if ( vscl > 0 )
        % Min of eps*vscale and the minimum non-zero coefficient:
        absCoeffs(~absCoeffs) = min( min(eps*vscl), ...
                                 min(absCoeffs(logical(absCoeffs(:)))) );
else
    % Add eps for zero CHEBTECHs:
    absCoeffs = absCoeffs + eps;
end

% Plot the coeffs:
scatter3(x, y, z, [], log10(absCoeffs(:)), 'filled')

% By default, set grid on
grid(gca, 'on')
colorbar
view(-30,35)

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h;
end

end