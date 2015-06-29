function cg = demo()
%DEMO   Return a random BVP CHEBGUI demo.
%   CG = CHEBGUI.DEMO() returns a random BVP CHEBGUI demo, on a CHEBGUI object
% format.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Find the folders which demos are stored in. The chebguiDemos folder lives in
% the trunk folder, find the path of the Chebfun trunk.
trunkPath = chebfunroot();

% Append directory information:
bvppath = fullfile(trunkPath, 'chebguiDemos', 'bvpdemos');

% Setup ODEs:
D = dir(bvppath);

% Find how many demos there are available:
numberOfDemos = length(D)-2; % First two entries are . and ..
    
% Create a random integer in 1:numberOfDemos. We avoid use of randi which was
% introduced in v2008b.
randnum = floor(numberOfDemos*rand) + 1;
selected = randnum + 2; % Need to shift the random integer by 2

% The full path to the demo file:
demoPath = fullfile(bvppath, D(selected).name);

% Convert the .guifile to a CHEBGUI object:
cg = chebgui.demo2chebgui(demoPath);

end
