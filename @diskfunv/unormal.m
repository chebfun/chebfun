function N = unormal()
%UNORMAL 

%   N = UNORMAL returns a DISKFUNV that is normal to the unit circle
%
%   See also  DOT, CURL


    x = diskfun(@(t,r) r.*cos(t)); 
    y = diskfun(@(t,r) r.*sin(t)); 

N = diskfunv(x,y);

end
