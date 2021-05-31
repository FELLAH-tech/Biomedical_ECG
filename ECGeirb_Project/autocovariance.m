function [g] = autocovariance(k,R_R2,mean)
    N = length(R_R2);
    sum = 0;
    for n=1:N-k
        sum = sum + (R_R2(n+k)-mean)*(R_R2(n)-mean);
    end
    g = sum/(N-k-1);
        
    
        