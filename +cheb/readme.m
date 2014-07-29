function readme
%README   Information about the CHEB package.

fn = mfilename('fullpath');
[pn,fn] = fileparts(fn);
type(fullfile(pn,'README.txt'))

end