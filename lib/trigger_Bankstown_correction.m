function header = trigger_Bankstown_correction(header)
%%%Bankstown trigger correction for delays

%%% send markers to workspace, then...
% 
% x=cell2mat(events(:,3));
% 
% diff_x=diff(x);
% figure();plot(diff_x);
% 
% [i,j]=find(diff_x>30);

events = squeeze(struct2cell(header.events))'; 

%ind=find(ismember(events(:,1),'D247'));
%ind=ind(1:2:end);
%events=events(ind+1,:);
ind=find([header.events{:,3}]==0);
ind=ind+1;
ind=ind';
header.events = header.events(ind,:);

