function out = sum(f, a, b)
%SUM   Definite integral of a CHEBFUN.
%   SUM(F) is the integral of a column CHEBFUN F over its domain of definition.
%
%   SUM(F, A, B) integrates a column CHEBFUN F over [A, B], which must be a
%   subdomain of F.domain. A and B may be scalars or CHEBFUN objects. In the
%   later case, S = SUM(F, A, B) returns a CHEBFUN S which satisfies
%                       B(s)
%                       /
%               S(s) =  | F(t) dt.
%                       /
%                     A(s)
%
%   SUM(F, DIM), where DIM is one of 1, 2, sums F over the dimension DIM. If F
%   is a column CHEBFUN and DIM = 1 or if F is a row CHEBFUN and DIM = 2 then
%   this integrates in the continuous dimension of F, as described above.
%   Otherwise, SUM(F, DIM) sums across the columns (rows) of the column (row)
%   CHEBFUN F.
%
% See also CUMSUM, DIFF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% An empty CHEBFUN has sum 0:
if ( isempty(f) )
    out = 0;
    return
end

% Integrate in the continuous dimension by default.
if ( f.isTransposed )
    dim = 2;
else
    dim = 1;
end

% Parse inputs:
doSubDomain = 0;
if ( nargin == 3 )
    doSubDomain = 1;
elseif ( (nargin == 2) && (numel(a) == 2) )
    % Support for sum(f, [a, b]):
    b = a(2);
    a = a(1);
    doSubDomain = 1;
elseif ( (nargin == 2) && (numel(a) == 1) )
    % Support for sum(f, dim):
    dim = a;
end

if ( xor(f.isTransposed, dim == 2) ) % Sum over the columns:
    % Call SUMCOLUMNS():
    out = sumColumns(f);
elseif ( doSubDomain )               % Integrate over a subdomain:
    % Call SUMSUBDOMAIN():
    out = sumSubDom(f, a, b);
else                                 % Integrate over the whole domain:
    % Call SUMFULLDOMAIN():
    out = sumFullDom(f);
end
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% SUM the columns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = sumColumns(f)
    % Sum the FUNS across the columns:
    for k = 1:numel(f.funs)
        f.funs{k} = sum(f.funs{k}, 2);
    end
    % Sum the impulses across the columns:
    f.impulses = sum(f.impulses, 2);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUM the whole domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sumFullDom(f)
    % Initialise out to be zero:
    out = 0;

    % Sum on each FUNs:
    for k = 1:numel(f.funs)
        out = out + sum(f.funs{k});
    end

    % Deal with impulses:
    if ( size(f.impulses, 3) > 1 )
        % Only add the delta functions (the integral of derivatives of delta
        % functions is always zero).
        out = out + sum(f.impulses(:,:,2));
    end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% SUM on a subdomain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sumSubDom(f, a, b)
    d1 = f.domain(1);
    d2 = f.domain(end);

    if ( isnumeric(a) && isnumeric(b) )

        % Validate the subdomain:
        if ( (a < d1) || (b > d2) )
            error('CHEBFUN:sum:ab', 'Not a valid subdomain.');
            
        elseif ( (a == d1) && (b == d2) )
            % Subdomain == original domain.
            out = sum(f);
            return
            
        end

        % It's faster to apply CUMSUM() than RESTRICT():
        out = cumsum(f);
        % Compute the area in between:
        out = feval(out, b) - feval(out, a);

    elseif ( isa(a, 'chebfun') && isa(b, 'chebfun') )

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the CHEBFUNs defining the limits:
        out = compose(b, out) - compose(a, out);

    elseif ( isa(b, 'chebfun') )

        % Validate the subdomain:
        if ( a < d1 )
            error('CHEBFUN:sum:a', 'Not a valid subdomain.');
        end

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the CHEBFUNs defining the limits:
        out = compose(b, out) - feval(out, a);

    elseif ( isa(a, 'chebfun') )

        % Validate the subdomain:
        if ( b > d2 )
            error('CHEBFUN:sum:b', 'Not a valid subdomain.');
        end

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the CHEBFUNs defining the limits:
        out = feval(out, b) - compose(a, out);
    end

end
