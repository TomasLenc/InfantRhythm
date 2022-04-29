function res=SNR(x, binmin, binmax)

res = zeros(size(x)); 
bl = zeros(size(x)); 

for d3=1:size(x,3)
    for i=1:length(x)
        dx1 = max(i-binmax,1); 
        dx2 = max(i-binmin,1); 
        dx3 = min(i+binmin,length(x)); 
        dx4 = min(i+binmax,length(x)); 

        bl(:,i,d3) = mean(x(:,[dx1:dx2,dx3:dx4],d3),2); 
    end
end

res = x-bl; 







