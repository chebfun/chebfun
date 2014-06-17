function out = sum(F, a, b)
%SUM   Definite integral of a CHEBFUN.
%   SUM(F) is the integral of a column CHEBFUN F over its domain of definition.
%
%   SUM(F, A, B), where A and B are scalars, integrates a column CHEBFUN F over
%   [A, B], which must be a subdomain of F.domain:
%
%                         B
%                         /
%               SUM(F) =  | F(t) dt.
%                         /
%                        A
%
%   SUM(F, A, B), where A and B are CHEBFUN objects, returns a CHEBFUN S which
%   satisfies
%
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Integrate in the continuous dimension by default.
if ( F(1).isTransposed )
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

if ( isempty(F) )
    % Empty chebfun has sum 0. (v4 compatibility).
    if ( xor(F(1).isTransposed, dim == 2) )
        out = F;
    else
        out = 0;
    end
    return
end

% Initialise the output:
out = cell(1, numel(F));

% Sum in the appropriate dimension:
if ( xor(F(1).isTransposed, dim == 2) ) % Sum over the columns:
    % Call SUMCOLUMNS():
    out = sumColumns(F);
elseif ( doSubDomain )                  % Integrate over a subdomain:
    % Loop over the columns:
    for k = 1:numel(F)
        % Call SUMSUBDOMAIN():
        out{k} = sumSubDom(F(k), a, b);
    end
    out = [out{:}];
else                                    % Integrate over the whole domain:
    % Loop over the columns:
    for k = 1:numel(F)
        % Call SUMFULLDOMAIN():
        out{k} = sumFullDom(F(k));
    end
    out = [out{:}];
end
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% SUM the columns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = sumColumns(f)
    if ( numel(f) == 1 )
        % Sum the FUNS across the columns:
        for k = 1:numel(f.funs)
            f.funs{k} = sum(f.funs{k}, 2);
        end
        % Sum the pointValues across the columns:
        f.pointValues = sum(f.pointValues, 2);
    else
        s = f(1);
        for k = 2:numel(f)
            s = s + f(k);
        end
        f = s;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUM the whole domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sumFullDom(f)
    % Initialise out to be zero:
    out = 0;

    % Sum on each FUNs:
    for k = 1:numel(f.funs)
        out = out + sum(f.funs{k});
    end
    
    % To avoid things like NaN + 1i*NaN:
    if ( isnan(out) )
        out = nan;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% SUM on a subdomain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sumSubDom(f, a, b)
    d1 = f.domain(1);
    d2 = f.domain(end);

    if ( isnumeric(a) && isnumeric(b) )

        % Validate the subdomain:
        if ( (a < d1) || (b > d2) )
            error('CHEBFUN:CHEBFUN:sum:sumSubDom:ab', 'Not a valid subdomain.');
            
        elseif ( (a == d1) && (b == d2) )
            % Subdomain == original domain.
            out = sum(f);
            return
            
        end

        % It's faster to apply CUMSUM() than RESTRICT():
        out = cumsum(f);
        % Compute the area in between:
        out = feval(out, b) - feval(out, a);
        return
    end
    
    % ComposeChebfuns doesn't allow scalar expansion, so we must loop over cols.
    numCols = numColumns(f);
    if ( numCols > 1 )
        out = chebfun();
        for k = numCols:-1:1
            fk = extractColumns(f,k);
            out(k) = sumSubDom(fk, a, b);
        end
        return
    end

    if ( isa(a, 'chebfun') && isa(b, 'chebfun') )

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the CHEBFUNs defining the limits:
        out = compose(b, out) - compose(a, out);

    elseif ( isa(b, 'chebfun') )

        % Validate the subdomain:
        if ( a < d1 )
            error('CHEBFUN:CHEBFUN:sum:sumSubDom:a', 'Not a valid subdomain.');
        end

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the CHEBFUNs defining the limits:
        out = compose(b, out) - feval(out, a);

    elseif ( isa(a, 'chebfun') )

        % Validate the subdomain:
        if ( b > d2 )
            error('CHEBFUN:CHEBFUN:sum:sumSubDom:b', 'Not a valid subdomain.');
        end

        % Compute the indefinite integral:
        out = cumsum(f);
        % Compose with the CHEBFUNs defining the limits:
        out = feval(out, b) - compose(a, out);
    end

end
