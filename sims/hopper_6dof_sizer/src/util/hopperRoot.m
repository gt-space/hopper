function root = hopperRoot(varargin)
%HOPPERROOT  Absolute path of the hopper_6dof_sizer project folder.
%   root = hopperRoot()                 project root
%   p    = hopperRoot('data', 'x.mat')  fullfile(root, 'data', 'x.mat')
%
%   Use this for files a script writes (data/, results/, outputs/) so they land
%   in the right folder whatever the current folder is. Files that are only
%   read can be found through the project path.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin > 0
    root = fullfile(root, varargin{:});
end
end
