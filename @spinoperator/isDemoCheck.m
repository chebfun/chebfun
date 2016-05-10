function out = isDemoCheck(char)
%ISDEMO   Check whether the string which was passed corresponds to a demo.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: This is not a good way to use exceptions handling. We should replace
% this with a more *cosmic* solution.
out = [];
try spinop(char);
catch 
    try spinop2(char);
    catch 
        try spinop3(char);
        catch 
            % It's not a demo because SPINOP/SPINOP2/SPINOP3 didn't recognise
            % it:
            out = 0;
        end
    end
end

% It's a demo:
if ( isempty(out) == 1 )
    out = 1;
end

end