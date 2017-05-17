function f = compose(f, op, g, data, pref)
%COMPOSE   Composition of BNDFUN objects.
%   H = COMPOSE(F, OP) returns a BNDFUN representing OP(F) where F is also a
%   BNDFUN object, and OP is a function handle.
%
%   H = COMPOSE(F, OP, G) returns OP(F, G) where F and G are BNDFUN objects, and
%   OP is a function handle. Note that the inputs are assumed to "make sense",
%   that is, F and G are assumed to have the same domain. No warning will be
%   given/error thrown if that is not the case, but the output of the method
%   will not be useful.
%
%   H = COMPOSE(F, G) returns a BNDFUN representing G(F), where both F and G are
%   also BNDFUN objects. If the range of F is not included in the domain of G,
%   an error is thrown. Notice that the domain of H will be the same as the
%   domain of F.
%
%   H = COMPOSE(F, OP, G, DATA, PREF) or H = COMPOSE(F, OP, [], DATA, PREF),
%   where F (and G) are BNDFUN objects, and OP is a function handle, uses the
%   constructor data in the structure DATA and the options passed by the
%   CHEBFUNPREF or preference structure PREF to build the returned BNDFUN. In
%   particular, if the underlying tech class supports it, one can use PREF to
%   alter the constructor's behavior to take advantage of the fact that F (and
%   possibly OP or G) has additional structure beyond just being an object.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin == 2 )
    if ( isa(op, 'bndfun') )  % Composition of two BNDFUN objects ( out = g(f) )
        
        % OP is a BNDFUN! Rename it to g for clarity
        g = op;
        
        % For composition to work with two BNDFUN objects, the range of F must
        % be within the domain of G. However, we allow this method to assume
        % that is the case.
        %
        % When we compose F and G, we need to find the appropriate part of the
        % domain of G on which the range of F lies. Since we then want to work
        % with the onefun of G, we apply the inverse map for the domain of G to
        % find the (sub)interval on [-1, 1] on which the range of F (and hence
        % its ONEFUN) lives.
        fMapped  = g.mapping.Inv(f.onefun);
        % Try to do the composition, will only work if the range of F lies in
        % the domain of G.
        f.onefun = compose(fMapped, g.onefun);
        
    else                     % Standard composition ( out = OP(F) )
        % OP is an operator!
        f.onefun = compose(f.onefun, op);
    end
    
elseif ( nargin == 3 )          % out = OP(F,G). No preferences passed.
    % Third argument passed must have been a bndfun. Compose the onefun of f
    % with the onefun of g:
    f.onefun = compose(f.onefun, op, g.onefun);
    
else                            % out = OP(F,G). Data (and maybe prefs) passed.
    % Get preferences.
    if ( nargin < 5 )
        pref = chebfunpref();
    else
        pref = chebfunpref(pref);
    end

    % Call compose:
    if ( isempty(g) )
        f.onefun = compose(f.onefun, op, g, data, pref);
    else
        f.onefun = compose(f.onefun, op, g.onefun, data, pref);
    end
end

end
