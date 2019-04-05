% THIS SCRIPT BRINGS ALL FIGURES TO THE TOP

all_figure_handles = findall(0, 'type', 'figure');
for iCount = length(all_figure_handles):-1:1
      figure(all_figure_handles(iCount))
end
clear all_figure_handles icount

