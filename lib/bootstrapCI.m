function [ci_high, ci_low] = bootstrapCI(x, n_boot) 
% calculate confidence intervals for an estimate of the mean 
% vals is a matrix with dims: [subjects]x[values]  

n_sub = size(x,1); 
n_vals = size(x,2); 

boot_idx = randsample([1:n_sub], n_sub*n_vals*n_boot, true); 

boot_idx = reshape(boot_idx, n_boot, n_vals, n_sub); 

boots = nan(n_boot, n_vals); 

for i=1:n_vals
    
    samples = x(:,i); 
    
    boots(:,i) = mean(samples(squeeze(boot_idx(:,i,:))), 2); 
    
end

ci_high = prctile(boots, 97.5, 1); 
ci_low = prctile(boots, 2.5, 1); 


% 
% figure
% n = 1000; 
% h = histogram(boots(:,n)); 
% hold on 
% plot([mean(x(:,n)),mean(x(:,n))], [0,max(h.Values)], '--k', 'linew',2); 
% plot([ci_high(:,n),ci_high(:,n)], [0,max(h.Values)], '--r', 'linew',2); 
% plot([ci_low(:,n),ci_low(:,n)], [0,max(h.Values)], '--r', 'linew',2); 
% 
% 
% figure
% i = [1:n_vals, n_vals:-1:1]; 
% ci = [ci_high, fliplr(ci_low)]; 
% fill(i, ci, 'r', 'facealpha',0.3, 'LineStyle','none')
