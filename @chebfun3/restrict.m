function f = restrict(f, dom)
% RESTRICT    Restrict the domain of a CHEBFUN3.
%
%   F = RESTRICT(F, DOM) returns a 
%   - CHEBFUN3 on the domain DOM that approximates F on that domain if DOM 
%     is a nondegenarate cuboid,
%   - CHEBFUN2 that approximates F on the domain if DOM is a plane, and
%   - CHEBFUN that approximates F on the domain if DOM is a line.
%   DOM should be a vector of length 6 specifying the underlying cuboid.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(dom, 'double') )
    if ( numel(dom) == 6 )                   % Restrict to DOM. 
        xlen = diff(dom(1:2));
        ylen = diff(dom(3:4));
        zlen = diff(dom(5:6));
        
        if ( ( xlen == 0 ) && ( ylen == 0) && ( zlen == 0 ) ) 
            % DOM is a point.
            f = feval(f, dom(1), dom(3), dom(5));
            
        elseif ( ( xlen ~= 0 ) && ( ylen == 0 ) && ( zlen == 0 ) )
            % DOM is a vertical line (corresponding to a column).
            cols = restrict(f.cols, dom(1:2));
            rows = feval(f.rows, dom(3));
            tubes = feval(f.tubes, dom(5));
            core = chebfun3.txm(chebfun3.txm(f.core, rows, 2), tubes, 3);
            f = chebfun(cols*core);
            
        elseif ( ( xlen == 0 ) && ( ylen ~= 0 ) && ( zlen == 0 ) )
            % DOM is a horizontal line (corresponding to a row).
            cols = feval(f.cols, dom(1));
            rows = restrict(f.rows, dom(3:4));
            tubes = feval(f.tubes, dom(5));
            core = chebfun3.txm(chebfun3.txm(f.core, cols, 1), tubes, 3);
            f = chebfun(rows*core.');
            
        elseif ( ( xlen == 0 ) && ( ylen == 0 ) && ( zlen ~= 0 ) )
            % DOM is an oblique line (corresponding to a tube).
            cols = feval(f.cols, dom(1));
            rows = feval(f.rows, dom(3));
            tubes = restrict(f.tubes, dom(5:6));
            core = squeeze(chebfun3.txm(chebfun3.txm(f.core, cols, 1), ...
                rows, 2));
            f = chebfun(tubes*core);
            
        elseif ( ( xlen == 0 ) && ( ylen ~= 0 ) && ( zlen ~= 0 ) )
            % DOM is a horizontal plane (corresponding to a horizontal slice).
            cols = feval(f.cols, dom(1));
            rows = restrict(f.rows, dom(3:4));
            tubes = restrict(f.tubes, dom(5:6));
            core = squeeze(chebfun3.txm(f.core, cols, 1));
            f = chebfun2(tubes*core.'*rows.');
            
        elseif ( ( xlen ~= 0 ) && ( ylen == 0 ) && ( zlen ~= 0 ) )
            % DOM is a lateral plane (corresponding to a lateral slice).
            cols = restrict(f.cols, dom(1:2));
            rows = feval(f.rows, dom(3));
            tubes = restrict(f.tubes, dom(5:6));
            core = squeeze(chebfun3.txm(f.core, rows, 2));
            f = chebfun2(tubes*core.'*cols.');
             
        elseif ( ( xlen ~= 0 ) && ( ylen ~= 0 ) && ( zlen == 0 ) )
            % DOM is a frontal plane (corresponding to a frontal slice).
            cols = restrict(f.cols, dom(1:2));
            rows = restrict(f.rows, dom(3:4));
            tubes = feval(f.tubes,dom(5));
            core = chebfun3.txm(f.core, tubes, 3);
            f = chebfun2(rows*core.'*cols.');
            
        else
            % DOM is not degenerate
            f.cols = restrict(f.cols, dom(1:2));
            f.rows = restrict(f.rows, dom(3:4));
            f.tubes = restrict(f.tubes, dom(5:6));
            f.domain = dom;
        end
    else
        error('CHEBFUN:CHEBFUN3:restrict:domain', 'Domain not determined.');
    end
        
else
    error('CHEBFUN:CHEBFUN3:restrict:domain', 'Unrecognizable domain.');
end

end
