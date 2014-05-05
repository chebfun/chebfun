function cg = loadexample(guifile,exampleNumber)

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Find the folders which demos are stored in. The chebguiDemos folder lives in
% the trunk folder, find the path of the Chebfun trunk.
trunkPath = fileparts(which('chebguiwindow'));

% Append directory information
bvppath = fullfile(trunkPath,'chebguiDemos','bvpdemos');

% Setup ODEs
D = dir(bvppath);

% exampleNumber == -1 denotes a random example is required
if exampleNumber == -1
    numberOfDemos = length(D)-2; % First two entries are . and ..
    % This can be used in V5.0
%     selected = randi(numberOfDemos)+2;
    
    % Avoid use of randi which was introduced in v2008b
    % Create a random integer in 1:numberOfDemos
    randnum = floor(numberOfDemos*rand)+1;
    selected = randnum + 2; % Need to shift the random integer by 2
else
    selected = exampleNumber + 2;
end

demoPath = fullfile(bvppath,D(selected).name);
cg = loaddemos(guifile,demoPath);
