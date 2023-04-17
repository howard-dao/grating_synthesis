function [] = save_all_figs( save_plots_path )
% Saves all the open figures

% get list of figure handles
figs = findobj( allchild(0), 'flat', 'Type', 'figure' );

% save each one
for ii = 1:length(figs)
    save_fig_multiformat( figs(ii), save_plots_path, get(figs(ii), 'Name' ), true );
end

<<<<<<< HEAD
end
=======
end

>>>>>>> b44a255
