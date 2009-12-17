function[] = plot_solution(handles, solution, limits);
% plot_solution -- Plots a DG solution
%
% plot_solution(handles, solution, limits)
%
%     'handles' is a length-K vector of graphics handles that point to the plot
%     of the solution on each element (there are K elements). Then 'solution' is
%     an N x K array of the solution values to be plotted. This function updates
%     all the plots using the values in solution. It imposes the axis scale
%     given by 'limits'.

K = length(handles);
for q = 1:K;
  set(handles(q), 'ydata', solution(:,q));
end
axis(limits);
drawnow;
