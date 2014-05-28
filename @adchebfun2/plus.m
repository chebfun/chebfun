function h = plus ( f, g )
%+	  Plus.
%
% F + G adds ADchebfun2s F and G, or a scalar to a chebfun2 if either F or G
% is a scalar.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if isempty(f) || isempty(g)  % check for empty ADchebfun2
   h = ADchebfun2;  % just return an empty ADchebfun2. 
   return; 
end

if ( ~isa(f,'adchebfun2') ) % First argument was not an ADchebfun2
    % Swap arguments.
    h = plus(g, f);
    return
end

if isa(g,'adchebfun2')     % ADCHEBFUN2 + ADCHEBFUN2
    h = f;
    
    % Add the chebfun2s
    h.chebfun2 = f.chebfun2 + g.chebfun2;
    
    % Add the der fields
    h.der = f.der + g.der;
    
elseif ( isa(g, 'chebfun2') || isa(g, 'double') ) % ADCHEBFUN2 + DOUBLE/SCALAR
    f.chebfun2 = f.chebfun2 + g;
else
    error('ADCHEBFUN2:plus:type','Cannot add these two objects together');
end



% if ( isa(f,'chebfun2') && isa(g,'double') )  
% %% chebfun2 + double
%     if isempty(f) % check for empty chebfun2.
%         return
%     end
%     h = f;
%     h.scl = f.scl;
%     h.fun2 = plus(f.fun2,g);
% elseif ( isa(f,'double') && isa(g,'chebfun2') )  
% %% double + chebfun2
%     if isempty(g) % check for empty chebfun2.
%         return
%     end
%     h = g;
%     h.scl = g.scl;
%     h.fun2 = plus(f,g.fun2);
% elseif ( isa(f,'chebfun2') && isa(g,'chebfun2') )  
% %% chebfun2 + chebfun2
%     if isempty(g) % check for empty chebfun2.
%         h = f; return
%     end
%     if isempty(f) % check for empty chebfun2.
%         h = g; return
%     end
%     h = f;
%     h.scl = f.scl;
%     h.fun2 = plus(f.fun2,g.fun2);
%     

% end

end