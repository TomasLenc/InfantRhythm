function letswaveNewPath = importLetswave(varargin)

addpath('..')
par = getParams();

if nargin==1 && varargin{1}==-1
    warning('off')
    rmpath(genpath(par.letswave_path)); 
    fprintf('\nremoving letswave from path...\n\n'); 
    warning('on')
    return
else
    addpath(genpath(par.letswave_path));
    fprintf('\nadding LW6 to path...\n\n'); 
end
