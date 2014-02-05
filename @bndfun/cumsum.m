function [f, rval] = cumsum(f, k, dim)
%CUMSUM   Indefinite integral of a BNDFUN.
%   CUMSUM(F) is the indefinite integral of the BNDFUN F on an interval [a,b],
%   with the constant of integration chosen so that F(a) = 0.
%
%   CUMSUM(F, K) will compute the Kth indefinite integral with the constant of
%   integration chosen so that each intermediate integral evaluates to 0 at x=a.
%   Thus CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, K, 2) will take the Kth cumulative sum over the columns F an
%   array-valued BNFUN.
%
% See also DIFF, SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%
% Trivial case of an empty BNDFUN:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 || isempty(k) )
    % Compute first indefinite intergral by default
    k = 1;
end

if ( nargin < 3 )
    % Assume dim = 1 by default
    dim = 1;
end

f = addBreaksForCumSum(f);
if ( iscell(f) )
    
    % If f is a cell of two BNDFUNs which are obtained by adding a new
    % breakpoint at the midpoint of the original domain due to non-trivial
    % exponents at both endpoints, we integrate each BNDFUN individually.
    
    for j = 1:2
        
        % Rescaling factor, (b-a)/2, to the kth power
        rescaleFactork = (.5*diff(f{j}.domain))^k;
        
        % Assign the ONEFUN of the output to be the output of the CUMSUM method 
        % of the ONEFUN of the input. If we called CUMSUM with third argument 
        % equal to 2 (i.e. % dim = 2), we only want to compute the cumlative sum
        % over columns, in which case, we should not rescale the result.
        if ( dim == 1 )
            f{j}.onefun = cumsum(f{j}.onefun, k, dim)*rescaleFactork;
        else
            f{j}.onefun = cumsum(f{j}.onefun, k, dim);
        end
    end
    % Get rval of the last bndfun:
    rval = get(f{2}, 'rval');
else
    % Rescaling factor, (b-a)/2, to the kth power
    rescaleFactork = (.5*diff(f.domain))^k;
    
    % Assign the ONEFUN of the output to be the output of the CUMSUM method of the
    % ONEFUN of the input. If we called CUMSUM with third argument equal to 2 (i.e.
    % dim = 2), we only want to compute the cumlative sum over columns, in which
    % case, we should not rescale the result.
    if ( dim == 1 )
        f.onefun = cumsum(f.onefun, k, dim)*rescaleFactork;
    else
        f.onefun = cumsum(f.onefun, k, dim);
    end
    
%     % Shift F up or down. This is useful at the chebfun level to concatenate the
%     % piece making the entire function as continuous as possible.
%     if ( ~any( issing(f) ) )
%         
%         % Grab the indices correspond to infinite shift:
%         ind = isinf(shift);
%         
%         % Zero the infinite shift:
%         shift( ind ) = 0;
%         
%         % Shift:
%         f = f + shift - get(f, 'lval');
%     end
    % Get rval of f:
    rval = get(f, 'rval');
    
end

end