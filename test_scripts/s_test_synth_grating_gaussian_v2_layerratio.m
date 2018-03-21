% authors: bohan
% 
% script for testing the new, NEW synthesis pipeline object


clear; close all;

% dependencies
addpath(['..' filesep 'main']);                                             % main code
addpath(['..' filesep '45RFSOI']);                                          % 45rf functions
% imec codes
addpath( 'C:\Users\bz\Google Drive\research\popovic group\projects\grating synthesis\code\IMEC_2018_02_07_tapeout\' ); 

% initial settings
disc                = 10;
units               = 'nm';
lambda              = 1200;
index_clad          = 1.448; %1.0;
domain              = [2500, 800];      % useful
optimal_angle       = 20;             % still useful
coupling_direction  = 'down';
data_dir            = 'C:\Users\bz\Google Drive\research\popovic group\projects\grating synthesis\data';
data_filename       = 'lol.mat';
data_notes          = 'meh';
data_mode           = 'new';
n_workers           = 1;

% handle to grating cell generation function
% IMEC cell
% etch_depth          = 70;
% h_makeGratingCell   = @( synth_obj, period, fill_top, fill_bot, offset_ratio ) ...
%                         f_makeGratingCell_IMEC( synth_obj, period, fill_top, fill_bot, offset_ratio, etch_depth );

% make object
Q = c_synthGrating( 'discretization',   disc,       ...
                    'units',            units,      ...
                    'lambda',           lambda,     ...
                    'background_index', index_clad, ...
                    'domain_size',      domain,     ...
                    'optimal_angle',    optimal_angle,      ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers, ...
                    'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
            );
%             'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...

        
% run synthesize gaussian grating
% note units are in 'nm'
% 
% angle   = 20; 
MFD     = 10000;
DEBUG   = false;
Q       = Q.synthesizeGaussianGrating(MFD, DEBUG);



% directivity vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, 10*log10(Q.directivities_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Directivity (dB) vs. fill factors');
savefig('dir_v_ff.fig');
saveas(gcf, 'dir_v_ff.png');

% directivity BEFORE sweeping periods vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, 10*log10(Q.dir_b4_period_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Directivity (dB) BEFORE PERIOD SWEEP vs. fill factors');
savefig('dir_b4_period_v_ff.fig');
saveas(gcf, 'dir_b4_period_v_ff.png');



% directivity vs. fill, saturated
dir_v_fill_sat                                  = 10*log10(Q.directivities_vs_fills);
sat_thresh                                      = 20;                                   % threshold, in dB
dir_v_fill_sat( dir_v_fill_sat < sat_thresh )   = sat_thresh;

figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, dir_v_fill_sat );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Directivity (dB) (saturated) vs. fill factors');
savefig([ 'dir_v_ff_sat_' num2str(sat_thresh) '.fig']);
saveas(gcf, [ 'dir_v_ff_sat_' num2str(sat_thresh) '.png']);


% % plot the way jelena did
% % with fill factor bottom vs. 'layer ratio'
% [FILL_BOT, FILL_TOP] = meshgrid( Q.fill_bots, Q.fill_tops );
% layer_ratio          = FILL_TOP./FILL_BOT;
% layer_ratio( isinf(layer_ratio) | isnan(layer_ratio) ) = 50;
% % layer_ratio = log10(layer_ratio);
% 
% figure;
% surf( layer_ratio, FILL_TOP, 10*log10(Q.directivities_vs_fills) );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('bottom layer ratio'); ylabel('top fill factor');
% title('Directivity (dB) vs. top fill factor and layer ratio');
% savefig('dir_v_ff_layer_ratio.fig');
% saveas(gcf, 'dir_v_ff_layer_ratio.png');


% % plot but block out all the places that jelena's code doesn't sweep
% dir_v_fill_jelena                       = 10*log10(Q.directivities_vs_fills);
% dir_v_fill_jelena( layer_ratio > 1.4 )  = -100;
% 
% figure;
% imagesc( Q.fill_bots, Q.fill_tops, dir_v_fill_jelena );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('bottom fill factor'); ylabel('top fill factor');
% title('Directivity (dB) (only datapoints that jelena sweeps) vs. fill factors');
% % savefig([ 'dir_v_ff_sat_' num2str(sat_thresh) '.fig']);
% % saveas(gcf, [ 'dir_v_ff_sat_' num2str(sat_thresh) '.png']);
% 
% 
% % plot but only show the "normal" curve regime
% dir_v_fill_jelena                       = 10*log10(Q.directivities_vs_fills);
% dir_v_fill_jelena( layer_ratio < 0.95 | layer_ratio > 1.4 ) = -100;
% 
% figure;
% imagesc( Q.fill_bots, Q.fill_tops, dir_v_fill_jelena );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('bottom fill factor'); ylabel('top fill factor');
% title('Directivity (dB) (only datapoints on normal curve) vs. fill factors');
% 
% 
% % plot but only show the "inverted" curve regime
% dir_v_fill_jelena                       = 10*log10(Q.directivities_vs_fills);
% dir_v_fill_jelena( layer_ratio > 1.05 | (layer_ratio + FILL_BOT) < 0.9 | (layer_ratio + FILL_BOT) > 1.1 )  = -100;
% 
% figure;
% imagesc( Q.fill_bots, Q.fill_tops, dir_v_fill_jelena );
% colorbar; set( gca, 'ydir', 'normal' );
% xlabel('bottom fill factor'); ylabel('top fill factor');
% title('Directivity (dB) (only datapoints on inverted curve) vs. fill factors');


% angles vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, Q.angles_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Angles (deg) vs. fill factors');
savefig('angle_v_ff.fig');
saveas(gcf, 'angle_v_ff.png');

% scattering strength alpha vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, real(Q.scatter_str_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Scattering strength (real) vs. fill factors');
savefig('scatter_str_v_ff.fig');
saveas(gcf, 'scatter_str_v_ff.png');

% period vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, Q.periods_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title(['Period (' Q.units.name ') vs. fill factors']);
savefig('period_v_ff.fig');
saveas(gcf, 'period_v_ff.png');

% offset vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, Q.offsets_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Offset vs. fill factors');
savefig('offsets_v_ff.fig');
saveas(gcf, 'offsets_v_ff.png');

% k vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, real(Q.k_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Real k vs. fill factors');
savefig('k_real_v_ff.fig');
saveas(gcf, 'k_real_v_ff.png');

figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, imag(Q.k_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Imag k vs. fill factors');
savefig('k_imag_v_ff.fig');
saveas(gcf, 'k_imag_v_ff.png');


% offset (jelena's definition) vs. fills
offset_jelena = Q.offsets_vs_fills + ...
                repmat(Q.fill_bots.', 1, length(Q.fill_top_bot_ratio) ) - ...
                repmat(Q.fill_bots.', 1, length(Q.fill_top_bot_ratio) ) .* repmat(Q.fill_top_bot_ratio, length(Q.fill_bots), 1);
offset_jelena = mod( offset_jelena, 1.0 );
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, offset_jelena );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Offset (jelena''s def) vs. fill factors');
savefig('offsets_jelena_v_ff.fig');
saveas(gcf, 'offsets_jelena_v_ff.png');


% plot of the final designs

% final synthesized fills
figure;
plot( 1:length(Q.bot_fill_synth), Q.bot_fill_synth, '-o' ); hold on;
plot( 1:length(Q.top_bot_fill_ratio_synth), Q.top_bot_fill_ratio_synth, '-o' );
xlabel('cell #');
legend('bottom fill ratio', 'top/bottom fill ratio');
title('Final synthesized fill factor');
makeFigureNice();
savefig('final_fill_ratios.fig');
saveas(gcf, 'final_fill_ratios.png');

% final synthesized fills
figure;
plot( 1:length(Q.bot_fill_synth), Q.bot_fill_synth.*Q.period_synth, '-o' ); hold on;
plot( 1:length(Q.top_bot_fill_ratio_synth), Q.top_bot_fill_ratio_synth.*Q.bot_fill_synth.*Q.period_synth, '-o' );
xlabel('cell #');
legend('bottom fill', 'top fill');
title('Final synthesized fills (nm)');
makeFigureNice();
savefig('final_fills.fig');
saveas(gcf, 'final_fills.png');

% final synthesized offset
figure;
plot( 1:length(Q.offset_synth), Q.offset_synth, '-o' );
xlabel('cell #');
legend('offset ratio');
title('Final synthesized offset ratio');
makeFigureNice();
savefig('final_offsets.fig');
saveas(gcf, 'final_offsets.png');

% final synthesized period
figure;
plot( 1:length(Q.period_synth), Q.period_synth, '-o' );
xlabel('cell #');
legend(['period, ' units]);
title('Final synthesized period');
makeFigureNice();
savefig('final_periods.fig');
saveas(gcf, 'final_periods.png');

% final synthesized angles
figure;
plot( 1:length(Q.angles_synth), Q.angles_synth, '-o' );
xlabel('cell #');
legend('angle, deg');
title('Final synthesized angles');
makeFigureNice();
savefig('final_angles.fig');
saveas(gcf, 'final_angles.png');

% final synthesized scattering strength
figure;
plot( 1:length(Q.scatter_str_synth), Q.scatter_str_synth, '-o' );
xlabel('cell #');
legend(['scattering strength (units 1/' units]);
title('Final synthesized scattering strength');
makeFigureNice();
savefig('final_scattering_str.fig');
saveas(gcf, 'final_scattering_str.png');

% final synthesized scattering strength normalized comparison
figure;
plot( 1:length(Q.scatter_str_synth), Q.des_scatter_norm, '-o' ); hold on;
plot( 1:length(Q.scatter_str_synth), Q.scatter_str_synth./abs(max(Q.scatter_str_synth(:))), '-o' );
xlabel('cell #');
legend('desired scatter str', 'normalized scatter str');
title('Final synthesized scattering strength, normalized comparison');
makeFigureNice();
savefig('final_scattering_str_norm.fig');
saveas(gcf, 'final_scattering_str_norm.png');

% final synthesized directivity
figure;
plot( 1:length(Q.dir_synth), 10*log10(Q.dir_synth), '-o' );
xlabel('cell #');
legend('Directivity, dB');
title('Final synthesized directivity');
makeFigureNice();
savefig('final_dir.fig');
saveas(gcf, 'final_dir.png');









        