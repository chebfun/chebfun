function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum on [-1,1].
%   VALS = MINANDMAX(F) returns a 2-vector VALS = [MIN(F); MAX(F)] with the
%   global minimum and maximum of the CHEBTECH F on [-1,1].  If F is a
%   vector-valued CHEBTECH, VALS is a 2-by-N matrix, where N is the number of
%   columns of F.  VALS(1, K) is the global minimum of the Kth column of F on
%   [-1, 1], and VALS(2, K) is the global maximum of the same.
%
%   [VALS, POS] = MINANDMAX(F) returns also the 2-vector POS where the minimum
%   and maximum of F occur.
%
%   If F is complex-valued the absolute values are taken to determine extrema
%   but the resulting values correspond to those of the original function. That
%   is, VALS = FEVAL(F, POS) where [~, POS] = MINANDMAX(ABS(F)).
%
% See also MIN, MAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isreal(f) )
    realf = real(f);
    imagf = imag(f);
    h = realf.*realf + imagf.*imagf;
    h = simplify(h);
    [~, pos] = minandmax(h);
    vals = feval(f, pos);
    vals = vals(:,1:size(pos,2)+1:end);
    return
end

% Compute derivative:
fp = diff(f);

% Make the Chebyshev grid (used in minandmaxColumn).
xpts = f.chebpts(length(f));

% Check if f has only one column (i.e., is not vector-valued):
sizef2 = size(f, 2);

if ( sizef2 == 1)
    % If f is a single column, then things are simple.
    [vals, pos] = minandmaxColumn(f, fp, xpts);
    
else
    % If f has multiple columns, we must loop over them.
    
    % Initialise output vectors:
    pos = zeros(2, sizef2);
    vals = zeros(2, sizef2);
    
    % Store coefficients and values:
    values = f.values;
    coeffs = f.coeffs;
    vscale = f.vscale;
    values2 = fp.values;
    coeffs2 = fp.coeffs;
    vscale2 = fp.vscale;
    
    % Loop over columns:
    for k = 1:sizef2    
        % Copy data to a single CHEBTECH column:
        f.values = values(:,k);
        f.coeffs = coeffs(:,k);
        f.vscale = vscale(k);
        fp.values = values2(:,k);
        fp.coeffs = coeffs2(:,k);
        fp.vscale = vscale2(k);

        % Find max of this column:
        [vals(:,k), pos(:,k)] = minandmaxColumn(f, fp, xpts);
    end
end

end

function [vals, pos] = minandmaxColumn(f, fp, xpts)

    % Initialise output
    pos = [ 0; 0 ];
    vals = [ 0; 0 ];

    % Compute turning points:
    r = roots(fp);
    r = [ -1; r; 1 ];
    v = feval(f, r);

    % min
    [vals(1), index] = min(v);
    pos(1) = r(index);

    % Take the minimum of the computed minimum and the function values:
    [vmin, vindex] = min([ vals(1); f.values ]);
    if ( vmin < vals(1) )
        vals(1) = vmin;
        pos(1) = xpts(vindex - 1);
    end

    % max
    [vals(2), index] = max(v);
    pos(2) = r(index);

    % Take the maximum of the computed maximum and the function values:
    [vmax, vindex] = min([ vals(2); f.values ]);
    if ( vmax > vals(2) )
        vals(2) = vmax;
        pos(2) = xpts(vindex - 1);
    end

end
