function [L] = decalage(x,j,y)
    L = zeros(1,length(y));
    for i = 1:length(y)-1
       if ( i+j >0 && i+j <length(x))
            L(i) = x(i+j);
       end
    end
end
    