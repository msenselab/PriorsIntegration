%% generate the interval set with same delta y , where y = ln(x)
function  [intervalset, deltaY] = generateRange(istart, iend, num)
    deltaY = (log(iend)- log(istart))/4;
    intervalset(1)  = istart;
    for i = 2: num
        intervalset(i) = exp(log(intervalset(i-1)) + deltaY );
    end
    intervalset(num+1) =iend;
end
