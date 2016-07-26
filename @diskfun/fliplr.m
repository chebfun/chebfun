function varargout +AD0- fliplr(varargin)
+ACU-FLIPLR   Flip/reverse a DISKFUN in the -direction.
+ACU-   G +AD0- FLIPLR( F ) returns a DISKFUN G with the same domain as F but reversed+ADs-
+ACU-   that is, G(x,y) +AD0- F(ax,y), where the domain is +AFs-a, b, c, d+AF0-.
+ACU-
+ACU- See also DISKFUN/FLIPUD.

+ACU- Copyright 2016 by The University of Oxford and The Chebfun Developers.
+ACU- See http://www.chebfun.org/ for Chebfun information.

+AFs-varargout+AHs-1:nargout+AH0AXQ- +AD0- fliplr+AEA-separableApprox(varargin+AHs-:+AH0-)+ADs-

end
