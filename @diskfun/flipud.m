function varargout +AD0- flipud(varargin)
+ACU-FLIPUD   Flip/reverse a DISKFUN in the latitude-direction.
+ACU-   G +AD0- FLIPUD(F) returns a DISKFUN G with the same domain as F but reversed+ADs-
+ACU-   that is, G(x,y) +AD0- F(x, cy), where the domain is +AFs-a, b, c, d+AF0-.
+ACU-
+ACU- See also DISKFUN/FLIPLR, DISKFUN/FLIPDIM. 

+ACU- Copyright 2016 by The University of Oxford and The Chebfun Developers.
+ACU- See http://www.chebfun.org/ for Chebfun information.

+AFs-varargout+AHs-1:nargout+AH0AXQ- +AD0- flipud+AEA-separableApprox(varargin+AHs-:+AH0-)+ADs-

end
