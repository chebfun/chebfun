function h = mean(f, g)
%MEAN    Average or mean value of a CHEBFUN.
%   MEAN(F) is the mean value of the CHEBFUN F in its continuous dimension.
%   MEAN(F, G) is the average CHEBFUN between CHEBFUN objects F and G.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: This needs to be reworked once we allow unbounded domains.

if ( nargin == 1 )

    if ( ~f.isTransposed )

        % Check to see if the domain is unbounded:
        infEnds = isinf(f.domain([1, end]));

        if ( ~any(infEnds) )
            % On a bounded domain, things are easy:
            h = sum(f)/diff(f.domain([1, end]));

        elseif ( all(infEnds == [1, 0]) )
            % If unbounded to the left, take the left impulse as the mean:
            h = f.impulses(1);

        elseif ( all(infEnds == [0, 1]) )
            % If unbounded to the right, take the right impulse as the mean:
            h = f.impulses(end);

        elseif ( all(infEnds == [1, 1]) )
            % If doubly unbounded...
            if ( abs(f.impulses(1) - f.impulses(end)) < get(f, 'vscale' ) )
                % Return the impulse at left/right if they are the same.
                h = f.impulses(1);
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
