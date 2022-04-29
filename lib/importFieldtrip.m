function importFieldtrip(varargin)

addpath('..')
par = getParams(); 

if nargin==1 & varargin{1}==-1
    warning('off')
    rmpath(genpath(par.fieldtrip_path)); 
else
    addpath(par.fieldtrip_path); 
    ft_defaults; 
end