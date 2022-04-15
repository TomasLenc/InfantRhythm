function [z_meterRel,z_meterUnrel,z,amps_meterRel,amps_meterUnrel,amps,frex_idx] = ...
    getZ(mX,frex,idxMeterRel,idxMeterUnrel,varargin)

% fix 1D data to be a row vector (remember: time must be the last dimension!)
if size(mX,1)>1 & size(mX,2)==1
    mX = mX'; 
end


if any(strcmpi(varargin,'fs')) & any(strcmpi(varargin,'N'))
    fs = varargin{find(strcmpi(varargin,'fs'))+1}; 
    N = varargin{find(strcmpi(varargin,'N'))+1}; 
    xstep = fs/N; 
elseif any(strcmpi(varargin,'xstep'))
    xstep = varargin{find(strcmpi(varargin,'xstep'))+1}; 
else
    error('you must to either specify both fs and N, or xstep in varargin!'); 
end
    
frex_idx = round(frex/xstep)+1; 

v = repmat({':'},ndims(mX),1); 
v{end} = frex_idx; 
amps = mX(v{:}); 

z = zscore(amps,[],ndims(mX)); 

v = repmat({':'},ndims(mX),1); 
v{end} = idxMeterRel; 
z_meterRel = mean(z(v{:}), ndims(mX)); 
amps_meterRel = mean(amps(v{:}), ndims(mX)); 

v = repmat({':'},ndims(mX),1); 
v{end} = idxMeterUnrel; 
z_meterUnrel = mean(z(v{:}), ndims(mX)); 
amps_meterUnrel = mean(amps(v{:}), ndims(mX)); 
