function varargout +AD0- flipdim(varargin)
+ACU-FLIPDIM   Flip/reverse a DISKFUN in a chosen direction.
+ACU-   G +AD0- FLIPDIM(F, DIM) returns a DISKFUN G with the same domain as F but
+ACU-   reversed in a direction, i.e., G(x,y)+AD0-F(x, cy). If DIM +AD0- 2 (default) then
+ACU-   G(x,y) +AD0- F(x, cy).  Otherwise DIM +AD0- 1 and G(x,y) +AD0- F(ax, y). The
+ACU-   domain of F is +AFs-a, b, c, d+AF0-.
+ACU- 
+ACU- See also DISKFUN/FLIPLR, DISKFUN/FLIPUD.

+ACU- Copyright 2016 by The University of Oxford and The Chebfun Developers.
+ACU- See http://www.chebfun.org/ for Chebfun information.

+AFs-varargout+AHs-1:nargout+AH0AXQ- +AD0- flipdim+AEA-separableApprox(varargin+AHs-:+AH0-)+ADs-

end
