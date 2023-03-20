function fsout = font(arg)
% FONT  Command window font size
%    font displays the font size
%    font + increases the font size
%    font - decreases the font size
%    font fs sets the font size
%    font(fs) sets the font size
%    fs = font gets the current font size

if nargin < 1 || isequal(arg,'+') || isequal(arg,'-')
   s = char(com.mathworks.services.FontPrefs.getCodeFont);
   if s(end-2) == '='
      fs = round(3/4*str2num(s(end-1)));
   else
      fs = round(3/4*str2num(s(end-2:end-1)));
   end
end
if nargin == 1
   if isequal(arg,'+')
      fs = fs+1;
   elseif isequal(arg,'-')
      fs = fs-1;
   elseif ischar(arg)
      fs = str2num(arg);
   else
      fs = arg;
   end
end
if nargout == 1
   fsout = fs;
else
   fprintf('font = %d\n',fs)
   com.mathworks.services.FontPrefs.setCodeFont( ...
      java.awt.Font('Lucida Sans Typewriter',java.awt.Font.BOLD,fs))
   com.mathworks.services.FontPrefs.setTextFont( ...
      java.awt.Font('SansSerif',java.awt.Font.PLAIN,fs))
end

