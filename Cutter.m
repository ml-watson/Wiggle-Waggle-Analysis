function y = Cutter(x,n,T)
% Type (t) --> changes type of cut
% T=1 --> Linear (ie n=2 every second)
% T=2 --> Log (ie n=2 takes 2nd, 4th, 8th, 16th)

if T == 1
    for i1 = (0:length(x)/n-1/n)
        y(i1+1) = x(1+i1*n);
    end
end

if T == 2
    % Check N value
    if n <= 1
        'Value of n must be greater than 1 for T=2'
        return
    end
    n1 = 0;
    i2 = 1;
    while (int64(n1 + n^(i2-1))) <= length(x)
       n1 = int64(n1 + n^(i2-1));
       y(i2) = x(n1);
       i2 = i2 + 1;
       if i2 >= log(length(x))/log(n)
           'Loop detected, Too many points'
           return
       end
    end
end
return