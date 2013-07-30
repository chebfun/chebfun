function out = sum(f, a, b)
%SUM    Definite integral of a chebfun.
%   SUM(F) is the integral of a chebfun F over the domain where it is defined.
%
%   SUM(F, A, B) integrates F over [A, B], which must be a subdomain of
%   F.domain. A and B may be scalars or chebfun objects. In the later case, S =
%   SUM(F, A, B) returns a chebfun S which satisfies
%                       B(s)
%                       /
%               S(s) =  | F(t) dt.
%                       /
%                     A(s)
%
% See also CUMSUM, DIFF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% An empty chebfun has sum 0:
if ( isempty(f) )
    out = 0;
    return
end

% Allow support for sum(f, [a, b]):
if ( nargin == 2 && numel(a) == 2)
    b = a(2);
    a = b(1);
end

if ( nargin == 1 )
    % Integrate over the whole domain:
    out = sumFullDom(f);
else
    % Integrate over a subdomain:
    out = sumSubDom(f, a, b);
end

end

%% SUM the whole domain:

function out = sumFullDom(f)
    % Initialise out to be zero:
    out = 0;

    % Sum on each funs:
    for k = 1:numel(f.funs)
        out = out + sum(f.funs{k});
    end

    % Deal with impulses:
    if ( size(f.impulses, 1) >= 2 )
        % Only add the delta functions (the integral of derivatives of delta
        % functions is always zero).
        out = out + sum(f.impulses(2,:));
    end

end

%% SUM a subdomain:

function out = sumSubDom(f, a, b)
    d1 = f.domain(1);
    d2 = f.domain(end);

    if ( isnumeric(a) && isnumeric(b) )

        % Validate the subdomain:
        if ( a < d1 || b > d2 )
            error('CHEBFUN:sum:ab', 'Not a valid subdomain.');
        elseif ( a == d1 && b == d2 )
            % Subdomain == original domain.
            out = sum(f);
            return
        end

        % It's faster to apply cumsum than to restrict:
        out = cumsum(f);
        % Compute the area in between:
        out = feval(out, b) - feval(out, a);

    elseif ( isa(a, 'chebfun') && isa(b, 'chebfun') )

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the chebfuns defining the limits:
        out = compose(out, b) - compose(out, a);

    elseif ( isa(b, 'chebfun') )

        % Validate the subdomain:
        if ( a < d1 )
            error('CHEBFUN:sum:a', 'Not a valid subdomain.');
        end

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the chebfuns defining the limits:
        out = compose(out, b) - feval(out, a);

    elseif ( isa(a, 'chebfun') )

        % Validate the subdomain:
        if ( b > d2 )
            error('CHEBFUN:sum:b', 'Not a valid subdomain.');
        end

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the chebfuns defining the limits:
        out = feval(out, b) - compose(out, a);
    end

end