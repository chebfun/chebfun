function h = mean(f, g)
%MEAN    Average or mean value of a CHEBFUN.
%   MEAN(F) is the mean value of the CHEBFUN F in its continuous dimension.
%   MEAN(F, G) is the average CHEBFUN between CHEBFUN objects F and G.
%
% See also SUM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )

    if ( ~f(1).isTransposed )

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

    else

        % Deal with the transposed case:
        h = transpose(mean(transpose(f)));
        % [TODO]: Is this right?

    end

else

    % Take the mean of two inputs:
    h = 0.5*(f + g);

end

end
