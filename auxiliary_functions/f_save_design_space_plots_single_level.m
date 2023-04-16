function [] = f_save_design_space_plots_single_level( synth_obj, fig_suffix, save_plots_path )
% Saves design space plots
%
% Inputs:
%   synth_obj
%       type: c_synthGrating object
%       desc: c_synthGrating object, after design space has been
%             generated
%   fig_suffix
%       type: string
%       desc: string to append to figure names, if you so desire. can be empty
%             string
%   save_plots_path
%       type: string
%       desc: path to save folder

% angles vs. fill
figure;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, ...
         synth_obj.sweep_variables.angles_vs_fill, '-o' );
xlabel('fill factor'); ylabel('angle');
title('Angles (deg) vs. fill factors');
figure_name = [ 'angle_v_ff' fig_suffix ];
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, figure_name );

% scattering strength alpha vs. fill
figure;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, ...
         real(synth_obj.sweep_variables.scatter_str_vs_fill), '-o' );
xlabel('fill factor'); ylabel('\alpha');
title('Scattering strength vs. fill factors');
figure_name = [ 'scatter_str_v_ff' fig_suffix ];
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, figure_name );

% predicted best MFD vs. fill for uniform designs
figure;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, ...
         2./(1.4625*real(synth_obj.sweep_variables.scatter_str_vs_fill)), '-o' );
xlabel('fill factor'); ylabel('MFD');
title('predicted best MFD vs. fill factors for unif design');
figure_name = [ 'best_MFD_v_ff' fig_suffix ];
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, figure_name );

% period vs. fill
figure;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, ...
         synth_obj.sweep_variables.periods_vs_fill, '-o' );
xlabel('fill factor'); ylabel('period');
title(['Period vs. fill factors']);
figure_name = [ 'period_v_ff' fig_suffix ];
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, figure_name );

% k vs. fill
figure;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, ...
         real(synth_obj.sweep_variables.k_vs_fill), '-o' );
xlabel('fill factor'); ylabel('k (real)');
title('Real k vs. fill factors');
figure_name = [ 'k_real_v_ff' fig_suffix ];
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, figure_name );

figure;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, ...
         imag(synth_obj.sweep_variables.k_vs_fill), '-o' );
xlabel('fill factor'); ylabel('k (imag)');
title('Imag k vs. fill factors');
figure_name = [ 'k_imag_v_ff' fig_suffix ];
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, figure_name );

end

