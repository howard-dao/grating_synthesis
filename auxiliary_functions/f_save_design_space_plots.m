function [] = f_save_design_space_plots( synth_obj, fig_suffix, save_plots_path )
% Saves design space plots
% Currently assumes sweep was top/bot
%
% Inputs:
%   synth_obj
%       type: c_synthTwoLeveLGrating object
%       desc: c_synthTwoLeveLGrating object, after design space has been
%             generated
%   fig_suffix
%       type: string
%       desc: string to append to figure names, if you so desire. can be empty
%             string
%   save_plots_path
%       type: string
%       desc: path to save folder

% % directivity vs. fill
% figure('name', [ 'dir_v_ff' fig_suffix ]);
% imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, 10*log10(synth_obj.sweep_variables.directivities_vs_fills) );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('top fill factor'); ylabel('bottom fill factor');
% title('Directionality (dB) vs. fill factors');
% % figure_name = [ 'dir_v_ff' fig_suffix ];
% % save_fig_multiformat( gcf, save_plots_path, figure_name );
% 
% % directivity BEFORE sweeping periods vs. fill
% try
%     figure('name', [ 'dir_b4_period_v_ff' fig_suffix ]);
%     imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, 10*log10(synth_obj.sweep_variables.dir_b4_period_vs_fills) );
%     colorbar; set( gca, 'ydir', 'normal' );
%     xlabel('top fill factor'); ylabel('bottom fill factor');
%     title('Directionality (dB) BEFORE PERIOD SWEEP vs. fill factors');
%     % figure_name = [ 'dir_b4_period_v_ff' fig_suffix ];
%     % save_fig_multiformat( gcf, save_plots_path, figure_name );
% catch
%     fprintf('no dir b4 period sweep to plot, skipping\n');
% end
% 
% % angles vs. fill
% figure('name', [ 'angle_v_ff' fig_suffix ]);
% imagesc( synth_obj.sweep_variables.fill_tops, ...
%          synth_obj.sweep_variables.fill_bots, ...
%          synth_obj.sweep_variables.angles_vs_fills );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('top fill factor'); ylabel('bottom fill factor');
% title('Angles (deg) vs. fill factors');
% % figure_name = [ 'angle_v_ff' fig_suffix ];
% % save_fig_multiformat( gcf, save_plots_path, figure_name );
% 
% % scattering strength (imaginary) alpha vs. fill
% figure('name', [ 'scatter_str_imag_v_ff' fig_suffix ]);
% imagesc( synth_obj.sweep_variables.fill_tops, ...
%          synth_obj.sweep_variables.fill_bots, ...
%          imag(synth_obj.sweep_variables.scatter_str_vs_fills) );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('top fill factor'); ylabel('bottom fill factor');
% title('Scattering strength (imaginary) vs. fill factors');
% % figure_name = [ 'scatter_str_imag_v_ff' fig_suffix ];
% % save_fig_multiformat( gcf, save_plots_path, figure_name );
% 
% % scattering strength alpha vs. fill
% figure('name', [ 'scatter_str_v_ff' fig_suffix ]);
% imagesc( synth_obj.sweep_variables.fill_tops, ...
%          synth_obj.sweep_variables.fill_bots, ...
%          real(synth_obj.sweep_variables.scatter_str_vs_fills) );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('top fill factor'); ylabel('bottom fill factor');
% title('Scattering strength (real) vs. fill factors');
% % figure_name = [ 'scatter_str_v_ff' fig_suffix ];
% % save_fig_multiformat( gcf, save_plots_path, figure_name );
% 
% % period vs. fill
% figure('name', [ 'period_v_ff' fig_suffix ]);
% imagesc( synth_obj.sweep_variables.fill_tops, ...
%          synth_obj.sweep_variables.fill_bots, ...
%          synth_obj.sweep_variables.periods_vs_fills );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('top fill factor'); ylabel('bottom fill factor');
% title(['Period vs. fill factors']);
% % figure_name = [ 'period_v_ff' fig_suffix ];
% % save_fig_multiformat( gcf, save_plots_path, figure_name );
% 
% % offset vs. fill
% figure('name', [ 'offsets_v_ff' fig_suffix ]);
% imagesc( synth_obj.sweep_variables.fill_tops, ...
%          synth_obj.sweep_variables.fill_bots, ...
%          synth_obj.sweep_variables.offsets_vs_fills );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('top fill factor'); ylabel('bottom fill factor');
% title('Offset vs. fill factors');
% % figure_name = [ 'offsets_v_ff' fig_suffix ];
% % save_fig_multiformat( gcf, save_plots_path, figure_name );
% 
% % k vs. fill
% figure('name', [ 'k_real_v_ff' fig_suffix ]);
% imagesc( synth_obj.sweep_variables.fill_tops, ...
%          synth_obj.sweep_variables.fill_bots, ...
%          real(synth_obj.sweep_variables.k_vs_fills) );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('top fill factor'); ylabel('bottom fill factor');
% title('Real k vs. fill factors');
% % figure_name = [ 'k_real_v_ff' fig_suffix ];
% % save_fig_multiformat( gcf, save_plots_path, figure_name );
% 
% figure('name', [ 'k_imag_v_ff' fig_suffix ]);
% imagesc( synth_obj.sweep_variables.fill_tops, ...
%          synth_obj.sweep_variables.fill_bots, ...
%          imag(synth_obj.sweep_variables.k_vs_fills) );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('top fill factor'); ylabel('bottom fill factor');
% title('Imag k vs. fill factors');
% % figure_name = [ 'k_imag_v_ff' fig_suffix ];
% % save_fig_multiformat( gcf, save_plots_path, figure_name );
% 
% 
% % plot altogether
% figsize = [ 600, 750 ];
% 
% % plot stuff i want
% figure('name', 'total_designspace', 'position', [400, 266, figsize ] );
% 
% % period vs fills
% s1 = subplot( 3, 2, 1 );
% 
% imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, 1e-3 * synth_obj.sweep_variables.periods_vs_fills );
% set(gca, 'ydir', 'normal'); axis image; colorbar;
% xlabel('D_u'); ylabel('D_l');
% title('\Lambda (\mum)');
% 
% % offsets vs fills
% s2 = subplot( 3, 2, 2 );
% 
% imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, 1e-3 *  synth_obj.sweep_variables.offsets_vs_fills );
% set(gca, 'ydir', 'normal'); axis image; colorbar;
% xlabel('D_u'); ylabel('D_l');
% title('O (\mum)');
% 
% % neff vs fills
% s3 = subplot( 3, 2, 3 );
% 
% imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, real(synth_obj.sweep_variables.k_vs_fills)./(2*pi/synth_obj.lambda) );
% set(gca, 'ydir', 'normal'); axis image; colorbar;
% xlabel('D_u'); ylabel('D_l');
% title('n_{eff}');
% 
% % alpha vs fills
% s4 = subplot( 3, 2, 4 );
% 
% imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, 1e3*imag(synth_obj.sweep_variables.k_vs_fills) );
% set(gca, 'ydir', 'normal'); axis image; colorbar;
% xlabel('D_u'); ylabel('D_l');
% title('\alpha (\mum^{-1})');
% 
% % theta vs fills
% s5 = subplot( 3, 2, 5 );
% 
% imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, synth_obj.sweep_variables.angles_vs_fills );
% set(gca, 'ydir', 'normal'); axis image; colorbar;
% xlabel('D_u'); ylabel('D_l');
% title(['\theta (' char(176) ')']);
% 
% % pdir vs fills
% s6 = subplot( 3, 2, 6 );
% 
% imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, 10*log10(synth_obj.sweep_variables.directivities_vs_fills) );
% set(gca, 'ydir', 'normal'); axis image; colorbar;
% xlabel('D_u'); ylabel('D_l');
% colormap(hot);
% title('10log_{10}(P_{dir})');
f_plot_design_space(synth_obj, fig_suffix);
save_all_figs(save_plots_path);

end

