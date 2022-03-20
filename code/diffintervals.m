function diffresult = diffintervals(intervals)
    diffresult = [ ];
    for i = 1 : length(intervals) - 1
        diffresult(i) = intervals(i+1) - intervals(i);
    end 
end 
