%% generate the interval set with same delta y , where y = ln(x)
function  intervalset = generateRangeDY(istart, deltaY, num)
    intervalset(1)  = istart;
    for i = 1: num
        intervalset(i+1) = exp( log(intervalset(i)) + deltaY );
    end
end
