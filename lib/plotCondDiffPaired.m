function plotCondDiffPaired(ax, values, col_lines, col_cond, varargin)
% plot lines connecting individual participant's tasks 

hold(ax,'on'); 

par = getParams(); 

mu = []; 
if any(strcmpi(varargin,'mu'))
    mu = varargin{find(strcmpi(varargin,'mu'))+1}; 
end

nCond = size(values,2); 

plot(ax, [1,2], values,'-',...
    'Color',col_lines, ...
    'LineWidth',1.5)

for iCond=1:nCond

    % plot zscore for each condition as point 
    plot(ax, iCond, values(:,iCond),'o',...
        'MarkerEdgeColor',col_cond(iCond,:),...
        'MarkerSize',5,...
        'LineWidth',1.5)    

    if ~isempty(mu) 
        % t-test against fixed value (e.g. cochlear model) 
        [H,P] = ttest(values(:,iCond), mu, 'tail', 'right'); 
        txt = text(ax, 1/4+(iCond-1)*2*1/4, 1.1, ...
                   num2str(round(P, 1, 'significant')), ...
                   'HorizontalAlignment','center', ...
                   'Units','Normalized'); 
    end
end

% paired t-test between conditions
if ~isempty(mu) 
    y = 1.3; 
else
    y = 1.1; 
end
[H,P] = ttest(values(:,1), values(:,2)); 
txt = text(ax, 0.5, y, ...
           ['p=', num2str(round(P, 1, 'significant'))], ...
           'HorizontalAlignment','center', ...
           'Units','Normalized'); 
txt.Color = par.col_neutral; 


ax.XLim = [0.5,2.5]; 
ax.Color = 'none'; 
ax.XTick = []; 
ax.XAxisLocation = 'origin'; 
ax.XTickLabel = []; 
ax.FontSize = par.fontsize; 
ax.LineWidth = 2;

