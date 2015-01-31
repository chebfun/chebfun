function f = cumsum(f, dim)
%CUMSUM   Indefinite integral of a BNDFUN.
%   CUMSUM(F) is the indefinite integral of the BNDFUN F on an interval [a,b],
%   with the constant of integration chosen so that F(a) = 0.
%
%   CUMSUM(F, 2) will take the cumulative sum over the columns F which is an
%   array-valued BNDFUN.
%
% See also DIFF, SUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Trivial case of an empty BNDFUN:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin < 2 )
    % Assume dim = 1 by default
    dim = 1;
end

if ( dim == 1 )
    
    % Integrate once to see if a cell is returned:
    g = cumsum(f.onefun);

    
    if ( iscell(g) )
        
        % If the integral of F.ONEFUN is a cell returned by CUMSUM@SINGFUN due 
        % to non-zero exponents at both endpoints, then we take care of each 
        % piece separately:
        
        dom = f.domain;
        emptyBndfun = bndfun();
        f = {emptyBndfun, emptyBndfun};
        
        % New domain and new map for the left piece:
        f{1}.domain = [dom(1) mean(dom)];
        f{1}.mapping = bndfun.createMap(f{1}.domain);
        
        % Rescaling for integration:
        f{1}.onefun = g{1}*diff(f{1}.domain);

        % New domain and new map for the right piece:
        f{2}.domain = [mean(dom) dom(2)];
        f{2}.mapping = bndfun.createMap(f{2}.domain);
        
        % Rescaling for integration:
        f{2}.onefun = g{2}*diff(f{2}.domain);
        
    else
        
        % Integrate once:
        f.onefun = (.5*diff(f.domain))*g;
        
    end
    
elseif ( dim == 2 )
    
    % When the third argument is 2, i.e. dim = 2, we compute the cumlative sum
    % over columns, in which case, no rescale is needed.
    
    f.onefun = cumsum(f.onefun, dim);
    
else
    error('CHEBFUN:BNDFUN:cumsum:input', ...
        'The third argument is unrecognizable.');
end

% Value of F at the right endpoint:
if ( iscell(f) )
    rVal = cellfun(@(f) get(f, 'rval'), f);
else
    rVal = get(f, 'rval');
end

end
