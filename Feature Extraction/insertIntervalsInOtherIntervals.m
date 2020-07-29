function new_intervals = insertIntervalsInOtherIntervals(old_intervals, ...
    intervals2add, plotFigure)

[row1, col1] = size(old_intervals);
[row2, col2] = size(intervals2add);

if col1~=2
    intervals2add = intervals2add';
end

if col2~=2
    old_intervals = old_intervals';
end

max_val = max(max([old_intervals; intervals2add]));
min_val = min(min([old_intervals; intervals2add]));


% get the logical vector for the old intervals
time_vector_intervals = 1:max_val;
log_vector_old_intervals = zeros(max_val,1);
for kk = 1:row1
    log_vector_old_intervals(old_intervals(kk,1):old_intervals(kk,2)) = 1;
end

% get the logical vector for the intervals to add
log_vector_intervals2add = zeros(max_val,1);
for kk = 1:row2
    log_vector_intervals2add(intervals2add(kk,1):intervals2add(kk,2)) = 1;
end

% get the logical vector for the new intervals
log_vector_new_intervals = double(log_vector_old_intervals & log_vector_intervals2add)';

if plotFigure
    figure(33)
    plot(time_vector_intervals,log_vector_old_intervals)
    hold on
    plot(time_vector_intervals,log_vector_intervals2add*2)
    plot(time_vector_intervals(min_val), log_vector_intervals2add(min_val)*2, '*k')
    plot(time_vector_intervals,log_vector_new_intervals*2.5)
    hold off
    ylim([0 3])
    xlim([min_val max_val])
    legend('old intervals', 'intervals to add', 'start intervals', 'new intervals')
end


% get the intervals for the new intervals

start_indexes = strfind(log_vector_new_intervals,[0 1])'+1;% add the zero for the cases in that the vector starts with zero
end_indexes = strfind(log_vector_new_intervals,[1 0])';
if log_vector_new_intervals(1)==1
    start_indexes = [1; start_indexes];
end
if log_vector_new_intervals(end)==1
    end_indexes = [end_indexes; length(log_vector_new_intervals)];
end


new_intervals = [start_indexes end_indexes];

end