function h = conv(f, g)
%CONV   Convolution of DELTAFUN objects.
%   H = CONV(F, G) produces the convolution of DELTAFUN objects F and G:
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = deltafun();
    return
end

if ( ~isa(f, 'deltafun') )
    % Then g must be a deltafun
    h = conv(g, f);
    return;
end

% f is definitely a deltafun at this point, extract deltafunctions of f:
f = simplifyDeltas(f);
g = simplifyDeltaf(g);
deltaMagF = f.deltaMag;
deltaLocF = f.deltaLoc;
deltaMagG = g.deltaMag;
deltaLocG = g.deltaLoc;

hFun = conv(f.funPart, g.funPart);

   if(isfimps)
        [m n] = size(fimps);
        % loop through the imps matrix
        for i = 1:m
            for j = 1:n
                if(abs(fimps(i,j)) > 100*eps)
                    % take appropriate derivative and shift the function
                    gshift = newdomain(diff(g,i-1),[g.ends(1)+f.ends(j) g.ends(end)+f.ends(j)]);
                    % pad with zero chebfuns on either side
                    l = chebfun( 0, [h.ends(1) gshift.ends(1) ] );
                    r = chebfun( 0, [gshift.ends(end) h.ends(end) ] );
                    gshift = chebfun( [ l; gshift; r ], [ h.ends(1) gshift.ends(1) gshift.ends(end) h.ends(end) ] );
                    % scale by the impulse value and add
                    h = h + fimps(i,j)*gshift;
                end                     
            end
        end
    end
    
    % if g has delta funtions, do the same as above 
    isgimps = any(any(abs(gimps)>100*eps));
    if(isgimps)
        [m n] = size(gimps);
        for i = 1:m
           for j = 1:n
               if(abs(gimps(i,j)) > 100*eps)
                   fshift = newdomain(diff(f,i-1),[f.ends(1)+g.ends(j) f.ends(end)+g.ends(j)]);
                   l = chebfun( 0, [h.ends(1) fshift.ends(1) ] );
                   r = chebfun( 0, [fshift.ends(end) h.ends(end) ] );
                   fshift = chebfun( [ l; fshift; r ], [ h.ends(1) fshift.ends(1) fshift.ends(end) h.ends(end) ] );
                   h = h + gimps(i,j)*fshift;
               end
           end
        end
    end
         
    % if both f and g have delta functions
    if(isfimps && isgimps)
        [m n] = size(fimps);
        [p q] = size(gimps);
        himps = zeros(m+p-1,length(h.ends));
        for i=1:m
            for j=1:n
                if(abs(fimps(i,j))>100*eps)
                    % find the locations of shifted ends in h.ends
                    [xx yy] = meshgrid(h.ends,f.ends(j)+g.ends);
                    idx = find(sum(~(abs(xx-yy)>100*eps)));
                    % place the scaled and shifted impulses in himps
                    himps(i:i+p-1,idx) = himps(i:i+p-1,idx) + fimps(i,j)*gimps;
                end
            end
        end
        % append delta functions to the imps of h
        h.imps = [ h.imps; himps ];    
    end
end


