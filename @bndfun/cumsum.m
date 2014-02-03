function f = cumsum(f, k, dim, shift)
%CUMSUM   Indefinite integral of a BNDFUN.
%   CUMSUM(F) is the indefinite integral of the BNDFUN F on an interval [a,b],
%   with the constant of integration chosen so that F(a) = 0.
%
%   CUMSUM(F, K) will compute the Kth indefinite integral with the constant of
%   integration chosen so that each intermediate integral evaluates to 0 at x=a.
%   Thus CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, K, 2) will take the Kth cumulative sum over the columns F an
%   array-valued BNDFUN.
%
%   CUMSUM(F, K, DIM, S) will shift F up by S. Note that this could be useful at
%   the CHEBFUN level to concatenate different pieces forming a countinuous
%   object.
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

if ( nargin < 4 )
    % Assume no need to shift:
    shift = 0;
end

% TODO: Please move this 'shift' operation back to the CHEBFUN level.

if ( dim == 1 )
    
    % Check if we need to add breakpoints:
    f = addBreaksForCumSum(f);
    
    if ( iscell(f) )
        
        % If f is a cell of two BNDFUNs which are obtained by adding a new
        % breakpoint at the midpoint of the original domain due to non-trivial
        % exponents at both endpoints, we integrate each BNDFUN individually.
        
        for j = 1:2
            
            % Rescaling factor, (b-a)/2, to the kth power
            rescaleFactork = (.5*diff(f{j}.domain))^k;
            
            % Assign the ONEFUN of the output to be the output of the CUMSUM method
            % of the ONEFUN of the input:
            
            f{j}.onefun = cumsum(f{j}.onefun, k, dim)*rescaleFactork;
        end
        
    else
        % Rescaling factor, (b-a)/2, to the kth power
        rescaleFactork = (.5*diff(f.domain))^k;
        
        % Assign the ONEFUN of the output to be the output of the CUMSUM method
        % of the ONEFUN of the input:
        f.onefun = cumsum(f.onefun, k, dim)*rescaleFactork;
        
        
        % Shift F up or down. This is useful at the chebfun level to concatenate the
        % piece making the entire function as continuous as possible.
        if ( ~any( issing(f) ) )
            
            % Grab the indice correspond to infinite shift:
            ind = isinf(shift);
            
            % Zero the infinite shift:
            shift( ind ) = 0;
            
            % Shift:
            f = f + shift - get(f, 'lval');
        end
        
    end
    
elseif ( dim == 2 )
    
    % When the third argument is 2, i.e. dim = 2, we compute the cumlative sum
    % over columns, in which case, no rescale is needed.
    
    f.onefun = cumsum(f.onefun, k, dim);
    
else
    error('CHEBFUN:BNDFUN:cumsum:input', ...
        'The third argument is unrecognizable.');
end

end