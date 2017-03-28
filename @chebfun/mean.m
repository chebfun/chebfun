function h = mean(f, g)
%MEAN    Average or mean value of a CHEBFUN.
%   MEAN(F) is the mean of the CHEBFUN F.
%   If F is a quasimatrix, the mean is along the continuous dimension.
%   MEAN(F, DIM), where DIM = 1 or 2, is the mean in dimension DIM.
%   MEAN(F, G), where G is a CHEBFUN, is the mean of F and G.
%
% See also SUM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( f(1).isTransposed )
    if nargin == 1
        g = 2;
    end

    if isscalar(g) && (g==1 || g==2)
        h = mean(f.', 3-g).';
    else
        h = mean(f.', g.').';
    end

else    % f is not transposed

    if (nargin == 1)

        % Check to see if the domain is unbounded:
        dom = domain(f);
        infEnds = isinf(dom([1, end]));

        if ( ~any(infEnds) )
            % On a bounded domain, things are easy:
            h = sum(f)/diff(dom([1, end]));

        elseif ( all(infEnds == [1, 0]) )
            % If unbounded to the left, take the left impulse as the mean:
            h = f.pointValues(1);

        elseif ( all(infEnds == [0, 1]) )
            % If unbounded to the right, take the right impulse as the mean:
            h = f.pointValues(end);

        elseif ( all(infEnds == [1, 1]) )
            % If doubly unbounded...
            if ( abs(f.pointValues(1) - f.pointValues(end)) < ...
                    max(vscale(f)*eps) )
                % Return the impulse at left/right if they are the same.
                h = f.pointValues(1);
            else
                % Or return a NaN if they are different:
                h = NaN;
            end
       end

    else  % case of 2 inputs

        if isscalar(g) && g == 2    
            h = sum(f,2)/size(f,2);

        elseif isscalar(g) && g == 1
            h = mean(f);

        else
        % Take the mean of two inputs:
            h = 0.5*(f + g);
        end
    end

end
