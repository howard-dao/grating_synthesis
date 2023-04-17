classdef c_synthTwoLevelGrating < c_synthGrating
% Synthesizes a high directivity 2-level grating in an arbitrary process
%
% Authors: bohan zhang
%
% this documentation needs to be updated
%
% Prerequisites/dependencies
%   - c_synthGrating.m
%   - c_twoLevelGratingCell.m
%   - the utility folder
%
%   The user should define their own custom grating unit cell
%   drawing function.
%   HOWEVER, this function MUST have the following inputs and outputs, IN
%   ORDER (this is outdated 12/30/22):
%       function GC = your_makeGratingCell_function( dxy, units, lambda, background_index, y_domain_size, ...
%                                          period, fill_top, fill_bot,
%                                          offset_ratio );
%           % makes and returns a c_twoLevelGratingCell object
%           inputs(outdated):
%               synth_obj
%               type: c_synthGrating object AS STRUCT
%               desc: c_synthGrating object AS STRUCT
%           period
%               type: double, scalar
%               desc: period of the grating cell
%           fill_top
%               type: double, scalar
%               desc: ratio of top layer to period
%           fill_bot
%               type: double, scalar
%               desc: ratio of bottom layer to period
%           offset_ratio
%               type: double, scalar
%               desc: ratio of bottom layer offset to period
%       outputs:
%           GC
%               type: c_twoLevelGratingCell object
%               desc: two level grating cell object
%
%   As an example:
%       see the default makeGratingCell() function that's included at the
%       end of this class definition. It makes a simple waveguide + air
%       cladding grating, no fancyness
%
%
%
% Inputs to constructor:
%   Inputs are name-value pairs:
%   'discretization'
%       type: double, scalar
%       desc: discretization along x and y, in units of 'units'
%
%   'units'
%       type: string
%       desc: name and scaling of spatial units, supports 'm'
%             (meters), 'mm' (millimeters), 'um' (microns), 'nm'
%             (nanometers)
%
%   'lambda'
%       type: double, scalar
%       desc: wavelength to solve at, in units 'units'
%
%   'background_index'
%       type: double, scalar
%       desc: value of background index
%
%   'domain_size'
%       type: 1x2 array, double
%       desc: domain size, [ y height, x length ]
%
%   'optimal_angle'
%       type: double, scalar
%       desc: desired output angle, in deg
%
%   'data_notes'
%       type: string
%       desc: optional verbose notes/descriptor for this simulation
%
%   'h_makeGratingCell'
%       type: function handle
%       desc: handle to grating drawing function
%             has to have these inputs in this order: 
%               dxy, background_index, y_domain_size,
%               period, fill_top, fill_bot, offset_ratio
%
% Examples:
%
%   % make synthesis object
%   synth_obj = c_synthTwoLevelGrating(   'discretization',    disc, ...
%                                       'units',             units,   ...
%                                       'lambda',            lambda, ...
%                                       'background_index',  background_index,    ...
%                                       'y_domain_size',     y_domain_size, ...
%                                       'optimal_angle',     optimal_angle, ...
%                                       'data_notes',        data_notes, ...
%                                       'coupling_direction', coupling_direction, ...
%                                       'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
%                                       );



    properties
        
        % inherits the following from c_synthGrating
        % discretization;         % dx and dy
        % units;                  % units, verbose, 'm' or 'mm', or 'um', or 'nm'
        %                         % has fields 'name' and 'scale'
        % lambda;                 % center wavelength
        % background_index;       % background index
        % y_domain_size;          % transverse domain size (y coordinate aka vertical dimension)
        % inputs;                 % saves input settings for user reference
        % start_time;             % time when object was created, 'YEAR-month-day hour-min-sec'
        % data_notes;             % verbose notes of what current sweep is doing     
        % h_makeGratingCell;      % handle to the grating cell making function
        % coupling_direction;     % either 'up' or 'down
        % optimal_angle;          % angle to optimize for, deviation from the normal, in deg.
        % sweep_variables;
        % synthesized_design;

        % fields saved in sweep_variables:
        % dimensions are top/bot ratio vs. bot fill
        % fill_tops
        % fill_bots
        % fill_top_bot_ratio
        % directivities_vs_fills
        % angles_vs_fills
        % scatter_str_vs_fills
        % periods_vs_fills
        % offsets_vs_fills - are these... ratios? or absolute values?
        % k_vs_fills
        % GC_vs_fills
        % dir_b4_period_vs_fills
                            
        % fields saved in synthesized_design:
        % dimensions are vs. cell #
        % dir                 
        % bot_fill             
        % top_bot_fill_ratio  
        % top_fill             
        % period               
        % offset        % absolute offset    
        % angles              
        % scatter_str         
        % k                 
        % GC               
        % des_scatter  
        % N
        % x_coords
        % y_coords
        % input_wg_type
        % use_min_feat_size     currently either true or false, doesn't save values
        
        % debug options
        % a struct that holds debugging options
        % currently saved fields:
        %   verbose
        debug_options;
        
        chosen_cells;
    end
    
    methods
        
        function obj = c_synthTwoLevelGrating(varargin)
            % Constructor
            % See top comments for input documentation
            
            % call c_synthGrating constructor
            obj = obj@c_synthGrating(varargin{:});
            
        end     % end constructor()

        

        % -----------------
        % Function generate_design_space()
        %
        %> Generates the design space vs. fill factors
        %> picking optimum offset and period for highest directivity and
        %> closest angle to desired
        %>
        %> This version runs versus top and bottom fill factors
        %>
        %> based on synthesizeGaussianGrating
        %>
        %> would be good to implement: option to save GC data or not
        %>
        %> Inputs:
        %>   fill_bots
        %>       type: double, array
        %>       desc: OPTIONAL Currently mostly for testing
        %>   fill_tops
        %>       type: double, array
        %>       desc: OPTIONAL Currently mostly for testing
        %   nworkers
        %       type: int, scalar
        %       desc: # of parallel workers to use (optional)
        %>   verbose
        %>       type: double, array
        %>       desc: OPTIONAL Currently mostly for testing, spits out
        %>             a bunch of stuff to the prompt
        function obj = generate_design_space( obj, fill_bots, fill_tops, nworkers, verbose )
            
            tic;
            fprintf('Sweeping fill factors for directivity and angle...\n');
            
            % set verbose options
            if ~exist('verbose', 'var')
                % default verbose off
                obj.debug_options.verbose = false;
            else
                obj.debug_options.verbose = verbose;
            end
            
            % set fill factors and offsets
            if ~exist('fill_bots', 'var')
                % default fill
                fill_bots   = fliplr( 0.025:0.025:0.975 );
            end
            if ~exist('fill_tops', 'var')
                % default fill
                fill_tops   = fliplr( 0.025:0.025:0.975 );
            end
            fill_top_bot_ratio  = [];
            guess_offset        = 0;
            
            % number of times to iterate over optimization loop
            n_optimize_loops = 2;
            
            % sort fills so they are in descending order
            fill_bots = sort( fill_bots, 'descend' );
            fill_tops = sort( fill_tops, 'descend' );
            
            % save fills and offsets
            obj.sweep_variables.fill_tops           = fill_tops;
            obj.sweep_variables.fill_bots           = fill_bots;
            obj.sweep_variables.fill_top_bot_ratio  = fill_top_bot_ratio;
            
            % initialize saving variables
            directivities_vs_fills  = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            angles_vs_fills         = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            periods_vs_fills        = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            offsets_vs_fills        = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill, this is offset ratio
            scatter_str_vs_fills    = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            k_vs_fills              = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            GC_vs_fills             = cell( length( fill_bots ), length( fill_tops ) );      % dimensions bot fill vs. top fill
            dir_b4_period_vs_fills  = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            prad_pin_vs_fills       = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            R_vs_fills              = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            
            % make grating cell, assuming both layers are filled
            waveguide = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               2*obj.discretization, ...
                                               1.0, 1.0, 0.0 );
            
            % run waveguide simulation
            % sim settings
            guess_n             = 1.0 * max( waveguide.N(:) );                                      % guess index. I wonder if there's a better guessk for this?
            guessk              = guess_n * 2*pi/obj.lambda;                                        % units rad/'units'
            num_wg_modes        = 5;
            BC                  = 0;                                                                % 0 = PEC
            pml_options_wg      = [0, 200, 20, 2];                                                  % now that I think about it... there's no reason for the user to set the pml options
            % run sim
            waveguide   = waveguide.runSimulation( num_wg_modes, BC, pml_options_wg, obj.k0, guessk );
            
            % update guessk (units rad/'units')
            guessk = waveguide.k;
            
            % grab waveguide k
            waveguide_k = waveguide.k;                                                              % units of rad/'units'
            
            % calculate analytical period which would approximately phase
            % match to desired output angle
            k0              = obj.background_index * ( 2*pi/obj.lambda );
            kx              = k0 * sin( (pi/180) * obj.optimal_angle );
            guess_period    = 2*pi/(waveguide_k - kx);                                              % units of 'units'
            
            % snap period to discretization
            guess_period = obj.discretization * round(guess_period/obj.discretization);
            
            % extend waveguide's ez for overlap
            [ waveguide, e_z_overlap_ext ]  = ...
                waveguide.stitch_E_field( waveguide.Phi, real(waveguide.k), round(guess_period/waveguide.domain_size(2)) );
            waveguide.E_z_for_overlap       = e_z_overlap_ext;
            
            % set grating solver settings
            num_modes   = 5;
            BC          = 0;                                                % 0 = PEC
            pml_options = [1, 100, 20, 2]; 
            OPTS        = struct( 'mode_to_overlap', e_z_overlap_ext );
            sim_opts    = struct('num_modes', num_modes, 'BC', BC, 'pml_options', pml_options);
            
            % initially start with waveguide GC
            guess_GC = waveguide;
                
            % first fill in the right side of the domain
            for i_ff_bot = 1:length( fill_bots )

                % print iteration
                fprintf('Right side of domain, fill iteration %i of %i\n', i_ff_bot, length(fill_bots) );

                % Optimize period and offset
                if fill_bots(i_ff_bot) < 1 && fill_tops(1) < 1
                    % Only run optimization if theres a perturbation of
                    % both layers
                    for ii = 1:n_optimize_loops
                        % run optimization loop
                        [ ~, ...
                          guess_period, ...
                          guess_offset, ...
                          directivities_vs_fills( i_ff_bot, 1 ), ...
                          angles_vs_fills( i_ff_bot, 1 ), ...
                          scatter_str_vs_fills( i_ff_bot, 1 ), ...
                          guess_GC, ...
                          guessk, ...
                          dir_b4_period_vs_fills( i_ff_bot, 1 ) ...
                          ] = ...
                            obj.optimize_period_offset( guess_offset, ...
                                                      fill_tops(1), ...
                                                      fill_bots(i_ff_bot), ...
                                                      guess_period,...
                                                      guessk, ...
                                                      sim_opts, ...
                                                      guess_GC );
                    end

                    % save the optimized period, offset, GC, and k
                    periods_vs_fills( i_ff_bot, 1 ) = guess_period;
                    offsets_vs_fills( i_ff_bot, 1 ) = guess_offset;
                    GC_vs_fills{ i_ff_bot, 1 }      = guess_GC;
                    k_vs_fills( i_ff_bot, 1 )       = guessk;
                    [~,R_vs_fills( i_ff_bot, 1 )]   = guess_GC.estimate_reflection(obj.background_index);
                    
                    if strcmp( obj.coupling_direction, 'up' )
                        prad_pin_vs_fills( i_ff_bot, 1 ) = guess_GC.P_rad_up/guess_GC.P_in;
                    else
                        prad_pin_vs_fills( i_ff_bot, 1 ) = guess_GC.P_rad_down/guess_GC.P_in;
                    end

                else
                    % at least one of the layers is not perturbed, so there is no optimization to run
                    % save dummy values
                    directivities_vs_fills( i_ff_bot, 1 )   = 1;
                    angles_vs_fills( i_ff_bot, 1 )          = 0;
                    scatter_str_vs_fills( i_ff_bot, 1 )     = 0;
                    periods_vs_fills( i_ff_bot, 1 )         = guess_period;
                    offsets_vs_fills( i_ff_bot, 1 )         = 0;
                    k_vs_fills( i_ff_bot, 1 )               = guessk;
                    GC_vs_fills{ i_ff_bot, 1 }              = waveguide;
                    dir_b4_period_vs_fills( i_ff_bot, 1 )   = 1;
                    prad_pin_vs_fills( i_ff_bot, 1 )        = 0;

                end     % end if fill_top < 1

                toc;

            end     % end initial domain sweep
            
            
            % for parallel processing, grab these variables before entering
            % parfor
            offsets_vs_fills_1 = offsets_vs_fills(:,1);
            periods_vs_fills_1 = periods_vs_fills(:,1);
            k_vs_fills_1       = k_vs_fills(:,1);
            GC_vs_fills_1      = GC_vs_fills(:,1);
            % calc number of loops, also necessary apparently for parfor
            n_fill_bots    = length(fill_bots);
            n_fill_tops    = length(fill_tops);

            % start up parallel pool
            if ~exist('nworkers', 'var')
                obj = obj.start_parpool();
            else
                obj = obj.start_parpool(nworkers);
            end
            
            
            % now fill in the rest of the domain
            parfor i_ff_bot = 1:n_fill_bots
                % For each bottom fill factor
    
                fprintf('Main parfor iteration %i of %i\n', i_ff_bot, n_fill_bots);
                
                % grab starting guess period and k
                guess_period    = periods_vs_fills_1( i_ff_bot );
                guessk          = k_vs_fills_1( i_ff_bot );
                guess_GC        = GC_vs_fills_1{ i_ff_bot };
                guess_offset    = offsets_vs_fills_1( i_ff_bot );

                % grab bottom fill
                fill_bot = fill_bots( i_ff_bot );

                for i_ff_top = 2:n_fill_tops
                    % for each top/bottom fill factor ratio

                    fprintf('Fill factor ratio %i of %i, main parfor iteration %i of %i\n', i_ff_top, n_fill_tops, i_ff_bot, n_fill_bots);
                        
                    fill_top = fill_tops(i_ff_top);
                    
                    % Optimize period and offset
%                     fill_top = fill_top_bot_ratio(i_ff_ratio) * fill_bot;
                    if fill_top < 1
                        % Only run optimization if theres a perturbation
                        
                        for ii = 1:n_optimize_loops
                            % run optimization loop
                            [ ~, ...
                              guess_period, ...
                              guess_offset, ...
                              best_directivity, ...
                              best_angle, ...
                              best_scatter_str, ...
                              guess_GC, ...
                              guessk, ...
                              dir_b4_period_vs_fill ...
                              ] = ...
                                obj.optimize_period_offset( guess_offset, ...
                                                          fill_top, ...
                                                          fill_bot, ...
                                                          guess_period,...
                                                          guessk, ...
                                                          sim_opts, ...
                                                          guess_GC );
                        end
                            
                        % save data
                        periods_vs_fills( i_ff_bot, i_ff_top )        = guess_period;
                        offsets_vs_fills( i_ff_bot, i_ff_top )        = guess_offset;
                        directivities_vs_fills( i_ff_bot, i_ff_top )  = best_directivity;
                        angles_vs_fills( i_ff_bot, i_ff_top )         = best_angle;
                        scatter_str_vs_fills( i_ff_bot, i_ff_top )    = best_scatter_str;
%                         GC_vs_fills{ i_ff_bot, i_ff_ratio }             = best_GC;
                        k_vs_fills( i_ff_bot, i_ff_top )              = guessk;
                        dir_b4_period_vs_fills( i_ff_bot, i_ff_top )  = dir_b4_period_vs_fill;
                        [~,R_vs_fills( i_ff_bot, i_ff_top )]          = guess_GC.estimate_reflection(obj.background_index);
                        if strcmp( obj.coupling_direction, 'up' )
                            prad_pin_vs_fills( i_ff_bot, i_ff_top ) = guess_GC.P_rad_up/guess_GC.P_in;
                        else
                            prad_pin_vs_fills( i_ff_bot, i_ff_top ) = guess_GC.P_rad_down/guess_GC.P_in;
                        end
                            
                    else
                        % at least one of the layers is not perturbed, so there is no optimization to run
                        % save dummy values
                        directivities_vs_fills( i_ff_bot, i_ff_top )    = 1;
                        angles_vs_fills( i_ff_bot, i_ff_top )           = 0;
                        scatter_str_vs_fills( i_ff_bot, i_ff_top )      = 0;
                        periods_vs_fills( i_ff_bot, i_ff_top )          = guess_period;
                        offsets_vs_fills( i_ff_bot, i_ff_top )          = 0;
                        k_vs_fills( i_ff_bot, i_ff_top )                = guessk;
%                         GC_vs_fills{ i_ff_bot, i_ff_ratio }               = waveguide;
                        dir_b4_period_vs_fills( i_ff_bot, i_ff_top )    = 1;
                        prad_pin_vs_fills( i_ff_bot, i_ff_top )         = 0;

                    end     % end if fill_top < 1
                end     % end for i_ff_top = ...         
            end     % end parfor i_ff_bot = ...
            
            % save variables to object
            obj.sweep_variables.directivities_vs_fills  = directivities_vs_fills;
            obj.sweep_variables.angles_vs_fills         = angles_vs_fills;
            obj.sweep_variables.scatter_str_vs_fills    = scatter_str_vs_fills;
            obj.sweep_variables.periods_vs_fills        = periods_vs_fills;
            obj.sweep_variables.offsets_vs_fills        = offsets_vs_fills;
            obj.sweep_variables.k_vs_fills              = k_vs_fills;
            obj.sweep_variables.dir_b4_period_vs_fills  = dir_b4_period_vs_fills;
            obj.sweep_variables.prad_pin_vs_fills       = prad_pin_vs_fills;
            obj.sweep_variables.R_vs_fills              = R_vs_fills;

            fprintf('Done generating design space\n');
            toc;
            
        end     % end generate_design_space()
        
        
        
        function obj = start_parpool( obj, nworkers )
            % convenience function for starting up a parallel pool
            
            % start a parallel pool session
            my_cluster = parcluster('local');                           % cores on compute node are "local"
            if getenv('ENVIRONMENT')                                    % true if this is a batch job
                my_cluster.JobStorageLocation = getenv('TMPDIR');       % points to TMPDIR
            end

            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if ~isempty(poolobj)
                % shut down previously made parallel pool
                delete(gcp('nocreate'));
            end
            % Automatically selects number of works to the max that the
            % cluster supports
            if nargin < 2
                parpool(my_cluster, my_cluster.NumWorkers);
            else
                parpool(my_cluster, nworkers); % temporary
            end
            
        end     % end start_parpool()

        
        
        function [  obj, best_period, best_offset, best_directivity, best_angle, best_scatter_str, ...
                    best_GC, best_k, dir_b4_period_vs_fill, DEBUG ] ...
                    = optimize_period_offset(obj, guess_offset, fill_top, fill_bot, guess_period, guessk, sim_opts, guess_gc )
            % for given fill ratios, optimizes period and offset for best angle/directivity
            %
            % Inputs:
            %   guess_offset
            %       type: scalar, double
            %       desc: guess offset position to start from, in units
            %             'units'
            %   fill_top
            %       type: scalar, double
            %       desc: ratio of top waveguide to period
            %   fill_bot
            %       type: scalar, double
            %       desc: ratio of bottom waveguide to period
            %   guess_period
            %       type: scalar, double
            %       desc: starting period to sweep from, in units 'units'
            %   guessk
            %       type: scalar, double
            %       desc: starting guess k
            %   sim_opts
            %       type: struct
            %       desc: structure with these fields:  
            %               num_modes
            %               BC
            %               pml_options
            %   guess_gc
            %       type: grating coupler object
            %       desc: initial grating coupler object to start with (for mode
            %       overlapping)
            %       
            %
            % Outputs:
            %   best_period
            %       type: scalar, double
            %       desc: optimal period
            %   best_offset
            %       type: scalar, double
            %       desc: absolute value of offset, in 'units'
            %   best_directivity
            %       type: scalar, double
            %       desc: optimal directivity
            %   best_angle
            %       type: scalar, double
            %       desc: optimal angle
            %   best_scatter_str
            %       type: scalar, double
            %       desc: optimal scatering strength
            %   best_GC
            %       type: c_twoLevelGratingCell
            %       desc: two level grating cell
            %   best_k
            %       type: scalar, double
            %       desc: optimal k
            %   dir_b4_period_vs_fill
            %       type: scalar, double
            %       desc: directivity b4 period sweep, mostly for debugging purposes 
            %   DEBUG
            %       type: struct
            %       desc: as name suggests, holds fields for debugging
            
            directivities = [];
            k_vs_offset   = [];
            angles        = [];
            GC_vs_offset  = {};

            % grab mode to overlap with
            OPTS = struct( 'mode_to_overlap', guess_gc.E_z_for_overlap );

            % new version of sweepign offsets, local search
            cur_offset   = guess_offset;
            offsets      = [];
            delta_offset = obj.discretization;  % start with positive delta offset
            i_offset     = 1;
            while true
               
                % make grating cell
                GC = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               guess_period, ...
                                               fill_top, ...
                                               fill_bot, ...
                                               cur_offset./guess_period );
                
                % run sim
                GC = GC.runSimulation( sim_opts.num_modes, sim_opts.BC, sim_opts.pml_options, obj.k0, guessk, OPTS );
                
                % save directivity
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    directivities( i_offset )   = GC.directivity;
                    angles( i_offset )          = GC.max_angle_up;
                else
                    % coupling direction is downwards
                    directivities( i_offset )   = 1./( GC.directivity );
                    angles( i_offset )          = GC.max_angle_down;
                end

                % update the guessk (units rad/'units')
                guessk                  = GC.k;
                k_vs_offset( i_offset ) = GC.k;
                
                % update the mode overlap field
                GC_vs_offset{i_offset}  = GC;
                OPTS.mode_to_overlap    = GC.E_z_for_overlap;
                
                offsets( i_offset ) = cur_offset;

                % Check if we need to switch directions
                if i_offset > 1
                    
                    if delta_offset > 0 && directivities( i_offset ) < directivities( i_offset-1 )
                        % Switch directions, we hit a max.
                        delta_offset            = -delta_offset;
                        OPTS.mode_to_overlap    = guess_gc.E_z_for_overlap;
                        cur_offset              = guess_offset;
                    elseif delta_offset < 0 && directivities( i_offset ) < max( directivities( 1:end-1) )
                        % Quit loop, we're done
                        break;
                    end
                    
                end
                
                % update for next loop
                cur_offset      = cur_offset + delta_offset;
                i_offset        = i_offset + 1;
                    
            end

            % DEBUG field
            DEBUG.GC_vs_offset = GC_vs_offset;

            % pick best offset, to feed into period sweep loop
            [ ~, indx_best_offset ]     = max( directivities );
%             best_offset_ratio           = offset_ratios( indx_best_offset );
            best_offset_ratio           = offsets( indx_best_offset )./guess_period;      % newer local search ver.
            best_offset_k               = k_vs_offset( indx_best_offset );
            best_offset_angle           = angles( indx_best_offset );
            
            % update mode overlap
            OPTS.mode_to_overlap    = GC_vs_offset{indx_best_offset}.E_z_for_overlap;

            % here's an output variable
            dir_b4_period_vs_fill = max( directivities );

            % now sweep periods
            % decide whether to sweep larger or smaller periods
            % based on the angle
            if best_offset_angle > obj.optimal_angle
                % only sweep smaller periods
                delta_period = -obj.discretization;
            else
                % only sweep larger periods
                delta_period = obj.discretization;
            end

            % init saving variables
            angles_vs_period    = []; 
            k_vs_period         = [];
            GC_vs_period        = {};
            periods             = [];
            
            % initial period sweep values
            guessk              = best_offset_k;
            periods(1)          = guess_period;
            GC_vs_period{1}     = GC_vs_offset{indx_best_offset};
            k_vs_period(1)      = best_offset_k;
            angles_vs_period(1) = best_offset_angle;
                       
            i_period    = 2;
            while true
             
                % verbose printing
                if obj.debug_options.verbose == true
                    fprintf('Sweeping period %i\n', i_period );
                end
                
                % update period
                guess_period = guess_period + delta_period;
                
                % make grating cell
                GC = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               guess_period, ...
                                               fill_top, ...
                                               fill_bot, ...
                                               best_offset_ratio );

                % run sim
                GC = GC.runSimulation( sim_opts.num_modes, sim_opts.BC, sim_opts.pml_options, obj.k0, guessk, OPTS );

                % save angle
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    angles_vs_period( i_period ) = GC.max_angle_up;
                else
                    % coupling direction is downwards
                    angles_vs_period( i_period ) = GC.max_angle_down;
                end

                % update for next iteration
                periods(i_period)       = guess_period;
                GC_vs_period{i_period}  = GC;
                k_vs_period(i_period)   = GC.k;
                guessk                  = GC.k;
                OPTS.mode_to_overlap    = GC.E_z_for_overlap;
                
                % check for exit condition (if error in angle gets worse)
                cur_angle_err   = abs( angles_vs_period( i_period ) - obj.optimal_angle );
                prev_angle_err  = abs( angles_vs_period( i_period-1 ) - obj.optimal_angle );
                if cur_angle_err > prev_angle_err
                    % optimization over, break
                    break;
                end
                
                i_period = i_period + 1;

            end     % end period sweep

            % DEBUG field
            DEBUG.GC_vs_period = GC_vs_period;
            
            % pick best period
            [angle_error, indx_best_period] = min( abs( obj.optimal_angle - angles_vs_period ) );
            best_period_k                   = k_vs_period( indx_best_period );
            
            % return data
            % best offset is already set
            best_GC = GC_vs_period{ indx_best_period };
            
            if strcmp( obj.coupling_direction, 'up' )
                % coupling direction is upwards
                best_directivity    = best_GC.directivity;
                best_angle          = best_GC.max_angle_up;
                best_scatter_str    = best_GC.alpha_up;
            else
                % coupling direction is downwards
                best_directivity    = 1./best_GC.directivity;
                best_angle          = best_GC.max_angle_down;
                best_scatter_str    = best_GC.alpha_down;
            end
            
            best_period = periods( indx_best_period );
            best_k      = best_GC.k;
            best_offset = round( best_offset_ratio * best_period / obj.discretization ) .* obj.discretization;  % snap to grid
                    
        end     % end function optimize_period_offset()
              
        
        
        function obj = generate_final_design_apodized( obj, desired_field, input_wg_type, enforce_min_feat_size_func, ...
                                                        fill_top_override, fill_bot_override, fill_vs_top_bot_both )
            % function for generating the final synthesized design
            % parameters
            %
            % uses the top and bottom fills
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter, in units 'units'
            %   input_wg_type
            %       type: string
            %       desc: 'bottom' for body only or 'full' for both layers
            %             or 'none' for no layers
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: OPTIONAL INPUT
            %             A function that user makes which enforces min.
            %             feat size
            %             See the example function at the bottom of this
            %             file
            %   fill_top_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_bot_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_vs_top_bot_both
            %       type: string
            %       desc: OPTIONAL: pick design vs top fill, bottom fill,
            %               or both
            %             either 'top', 'bottom', or 'both'
            %             default to 'top'
            %   
            %
            % The final design is saved as:
            %     obj.synthesized_design           
            
            % calcualte desired scattering
            [ obj, xvec, alpha_des, desired_field ] = obj.calculate_desired_scattering( desired_field );
            
            % select cells to use
            obj = obj.choose_unit_cells( input_wg_type, enforce_min_feat_size_func, ...
                                        fill_top_override, fill_bot_override, fill_vs_top_bot_both );

                     
            % optimize start alpha
            [ obj, best_alpha_power ] = obj.optimize_start_alpha( xvec, alpha_des, obj.chosen_cells.chosen_dirs, ...
                                              obj.chosen_cells.chosen_bot_fills, obj.chosen_cells.chosen_top_fills, ...
                                              obj.chosen_cells.chosen_offsets, obj.chosen_cells.chosen_periods, ...
                                              obj.chosen_cells.chosen_angles, obj.chosen_cells.chosen_scatter_str, ...
                                              obj.chosen_cells.chosen_ks, desired_field );
                                          
            % DEBUG PLOTTING PREDICTED FIELD SHAPE
            [ ~, field_shape_prediction ] = obj.predict_overlap_for_optimization( ...
                                              xvec, alpha_des, obj.chosen_cells.chosen_dirs, ...
                                              obj.chosen_cells.chosen_bot_fills, obj.chosen_cells.chosen_top_fills, ...
                                              obj.chosen_cells.chosen_offsets, obj.chosen_cells.chosen_periods, ...
                                              obj.chosen_cells.chosen_angles, obj.chosen_cells.chosen_scatter_str, ...
                                              obj.chosen_cells.chosen_ks, desired_field, best_alpha_power );
%             figure('Name','predicted_field');
%             plot( xvec, alpha_des ); hold on;
%             plot( xvec, abs(desired_field).*max(alpha_des)./max(desired_field) );
%             plot( xvec, field_shape_prediction.*max(alpha_des)./max(field_shape_prediction) );
%             xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
%             legend('desired scattering strength', 'fiber mode', 'predicted field shape');
%             title('DEBUG predicted field shape');
%             makeFigureNice();  
                                          
            % synthesize final design
            [ obj, synthesized_des ] = obj.pick_final_datapoints( xvec, alpha_des, ...
                                 obj.chosen_cells.chosen_dirs, ...
                                 obj.chosen_cells.chosen_bot_fills, ...
                                 obj.chosen_cells.chosen_top_fills, ...
                                 obj.chosen_cells.chosen_offsets, ...
                                 obj.chosen_cells.chosen_periods, ...
                                 obj.chosen_cells.chosen_angles, ...
                                 obj.chosen_cells.chosen_scatter_str, ...
                                 obj.chosen_cells.chosen_ks, ...
                                 10.^best_alpha_power );
            obj.synthesized_design = catstruct( obj.synthesized_design, synthesized_des );  % in aux function
            % save input waveguide type
            obj.synthesized_design.input_wg_type = input_wg_type; 
            
            % build final index distribution
            obj = obj.build_final_index();
            
        end     % end function generate_final_design_apodized()
        
        function obj = generate_final_design_apodized_gaussian( obj, MFD, input_wg_type, ...
                enforce_min_feat_size_func, fill_top_override, fill_bot_override, fill_vs_top_bot_both )
            % Synthesizes an apodized grating that radiates desired Gaussian field profile
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter of gaussian, defined as 1/e^2 width of intensity
            %   input_wg_type
            %       type: string
            %       desc: 'full' or 'bottom' are currently supported
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: user defined function for enforcing minimum feature sizes
            %               currently this function must take 2 args - ( period, fill )   
            %   fill_top_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_bot_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_vs_top_bot_both
            %       type: string
            %       desc: OPTIONAL: pick design vs top fill, bottom fill,
            %               or both
            %             either 'top', 'bottom', or 'both'
            %             default to 'top'
            
            % generate a fiber gaussian mode
            [ obj, field_profile ] = obj.make_gaussian_profile( MFD );
            
            % generate final design, apodized
            if nargin < 5
                obj = obj.generate_final_design_apodized( field_profile, input_wg_type, enforce_min_feat_size_func );
            elseif nargin < 7
                obj = obj.generate_final_design_apodized( field_profile, input_wg_type, enforce_min_feat_size_func, ...
                                                            fill_top_override, fill_bot_override );
            else
                obj = obj.generate_final_design_apodized( field_profile, input_wg_type, enforce_min_feat_size_func, ...
                                                            fill_top_override, fill_bot_override, fill_vs_top_bot_both );
            end
            
            % save MFD
            obj.synthesized_design.MFD = MFD;
            
        end % end function generate_final_design_apodized_gaussian()
             
        function obj = generate_final_design_uniform( obj, desired_field, input_wg_type, enforce_min_feat_size_func, filltop, fillbot )
            % Synthesizes a uniform grating that has amplitude matching field_profile
            %
            % reworking this code?
            % inputs:
            %   desired_field
            %       type: double, array
            %       desc: field profile to amplitude match to. dimensions
            %             are amplitude vs. x coordinates
            %   input_wg_type
            %       type: string
            %       desc: 'full' or 'none' are currently supported
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: OPTIONAL, user defined function for enforcing minimum feature sizes
            %               currently this function must take 4 args - ( period, top_fill, bot_fill, offset )
            %   filltop
            %       type: double, scalar
            %       desc: OPTIONAL fill to use
            %   fillbot
            %       type: double, scalar
            %       desc: OPTIONAL fill to use
            
                    
            % pick the best fills
            % select subset of domain to use
            % dimensions are bot vs. top
            if nargin < 4
                [ obj, top_fills, bot_fills, directivities, ...
                           angles, periods, offsets, scatter_strs, ks ] ...
                           = obj.pick_design_quadrant( input_wg_type );
            else
                % override fills
                [ obj, top_fills, bot_fills, directivities, ...
                           angles, periods, offsets, scatter_strs, ks ] ...
                           = obj.pick_design_quadrant( input_wg_type , filltop, fillbot );
            end
                   
            
            % init saving vars
            % dims are bot fill vs top fill
            obj.sweep_variables.overlap_predict_vs_fill = zeros( size( top_fills ) );

            % calculate transmission, bot fill vs top fill
            t = directivities./(directivities+1);

            % for each fill
            for i_bot = 1:size( bot_fills, 1 )

                for i_top = 1:size( top_fills, 2 )

                    % grab current alpha
                    cur_alpha = scatter_strs( i_bot, i_top );

                    % calc length of grating to scatter 99% of light (1/e^4 power)
                    grat_len = 2/cur_alpha;     % in 'units'

                    % gen x coords
                    xvec = 0 : obj.discretization : grat_len - obj.discretization;

                    % integrate the picked scatter strengths to get the predicted
                    % field shape
                    scatter_str_vs_x                = cur_alpha .* ones( size(xvec) );
                    [obj, field_shape_prediction]   = obj.predict_field_shape( xvec, scatter_str_vs_x );

                    % calculate overlap with desired field (simple xcorr)
                    [ ~, overlap ]  = obj.xcorr_normalized( desired_field, field_shape_prediction );
                    obj.sweep_variables.overlap_predict_vs_fill(i_bot, i_top) = max( abs(overlap) );

                    % check for min feat size, if doesn't satisfy, set overlap to 0
                    if exist( 'enforce_min_feat_size_func', 'var' )
                        if ~enforce_min_feat_size_func( periods( i_bot, i_top ), ...
                                top_fills( i_bot, i_top ), ...
                                bot_fills( i_bot, i_top ), ...
                                offsets( i_bot, i_top ) )
                            obj.sweep_variables.overlap_predict_vs_fill(i_bot, i_top) = 0;
                        end
                    end

                end     % end for i_top

            end     % end for i_bot

            % calculate coupling
            obj.sweep_variables.coupling_predict_vs_fill = obj.sweep_variables.overlap_predict_vs_fill.*t;
            obj.sweep_variables.fill_top_predict = top_fills( 1, : );
            obj.sweep_variables.fill_bot_predict = bot_fills( :, 1 );

            % now pick the best datapoint
            [ ~, i_bot_best, i_top_best ] = f_get_max_matrix_elem( obj.sweep_variables.coupling_predict_vs_fill );

            % generate final design
            grat_len            = 4./scatter_strs( i_bot_best, i_top_best );
            xvec                = 0 : obj.discretization : grat_len - obj.discretization;
            alpha_des           = scatter_strs( i_bot_best, i_top_best ) .* ones( size( xvec ) );
            [ obj, synthesized_des ] = obj.pick_final_datapoints( xvec, alpha_des, ...
                                             directivities( i_bot_best, i_top_best ), ...
                                             bot_fills( i_bot_best, i_top_best ), ...
                                             top_fills( i_bot_best, i_top_best ), ...
                                             offsets( i_bot_best, i_top_best ), ...
                                             periods( i_bot_best, i_top_best ), ...
                                             angles( i_bot_best, i_top_best ), ...
                                             scatter_strs( i_bot_best, i_top_best ), ...
                                             ks( i_bot_best, i_top_best ) );
            obj.synthesized_design = synthesized_des;
            obj.synthesized_design.input_wg_type = input_wg_type;
            
            % build final index distribution
            obj = obj.build_final_index();
            
        end     % end generate_final_design_uniform()
          
        function obj = generate_final_design_uniform_selectedfills( obj, input_wg_type, filltop, fillbot )
            % generate a final uniform design from user selected fills
            % inputs:
            %   input_wg_type
            %       type: string
            %       desc: 'bottom', 'top', 'full' or 'none'. mainly for
            %       documentation purposes
            %   filltop
            %       type: double, scalar
            %       desc: top layer fill (duty cycle) to use
            %   fillbot
            %       type: double, scalar
            %       desc: bottom layer fill (duty cycle) to use
                               
            directivities   = interp2( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots, obj.sweep_variables.directivities_vs_fills, filltop, fillbot );
            angles          = interp2( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots, obj.sweep_variables.angles_vs_fills, filltop, fillbot ); 
            periods         = interp2( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots, obj.sweep_variables.periods_vs_fills, filltop, fillbot );
            offsets         = interp2( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots, obj.sweep_variables.offsets_vs_fills, filltop, fillbot );
            scatter_strs    = interp2( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots, obj.sweep_variables.scatter_str_vs_fills, filltop, fillbot );
            ks              = interp2( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots, obj.sweep_variables.k_vs_fills, filltop, fillbot );
            
            % generate final design
            grat_len            = 4./scatter_strs;
            xvec                = 0 : obj.discretization : grat_len - obj.discretization;
            alpha_des           = scatter_strs .* ones( size( xvec ) );
            [ obj, synthesized_des ] = obj.pick_final_datapoints( xvec, alpha_des, ...
                                             directivities, ...
                                             fillbot, ...
                                             filltop, ...
                                             offsets, ...
                                             periods, ...
                                             angles, ...
                                             scatter_strs, ...
                                             ks );
            obj.synthesized_design = synthesized_des;
            obj.synthesized_design.input_wg_type = input_wg_type;
            
            % build final index distribution
            obj = obj.build_final_index();

        end
    

        function obj = generate_final_design_uniform_gaussian( obj, MFD, input_wg_type, enforce_min_feat_size_func, filltop, fillbot )
            % Synthesizes a uniform grating that radiates desired Gaussian field profile
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter of gaussian, defined as 1/e^2 width of intensity
            %   input_wg_type
            %       type: string
            %       desc: 'full' or 'none' are currently supported
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: user defined function for enforcing minimum feature sizes
            %               currently this function must take 2 args - ( period, fill ) 
            %   filltop
            %       type: double, scalar
            %       desc: OPTIONAL fill threshold to use
            %   fillbot
            %       type: double, scalar
            %       desc: OPTIONAL fill threshold to use
            
            % generate a fiber gaussian mode
            [ obj, field_profile ] = obj.make_gaussian_profile( MFD );
            
            % generate final design, uniform
            if nargin > 4
                obj = obj.generate_final_design_uniform( field_profile, input_wg_type, enforce_min_feat_size_func, filltop, fillbot );
            else
                obj = obj.generate_final_design_uniform( field_profile, input_wg_type, enforce_min_feat_size_func );
            end
            
            % save MFD
            obj.synthesized_design.MFD = MFD;
        
        end     % end generate_final_design_uniform_gaussian()
        
        function obj = choose_unit_cells( obj, input_wg_type, enforce_min_feat_size_func, ...
                                        fill_top_override, fill_bot_override, fill_vs_top_bot_both, ...
                                        interp_on, plot_on )
            % added a new input - interp_on. optionally interpolates the
            % final chosen cells. default is false
            % plot_on - optional. decide whether to plot the chosen cells
            % or not. default is true
                                    
            % select subset of domain to use
            if nargin < 4
             [ obj, top_fills, bot_fills, directivities, ...
                       angles, periods, offsets, scatter_strs, ks ] ...
                       = obj.pick_design_quadrant( input_wg_type );
            else
             % override the fills as user desired
             [ obj, top_fills, bot_fills, directivities, ...
                       angles, periods, offsets, scatter_strs, ks ] ...
                       = obj.pick_design_quadrant( input_wg_type, fill_top_override, fill_bot_override );
            end
            
            if nargin < 7
                interp_on = false;
            end
            
            if nargin < 8
                plot_on = true;
            end
            
            % default to using top fill
            if nargin < 6
                fill_vs_top_bot_both = 'top';
            end
            
            % remove unit cells that don't comply with min feature rules
            for i_bot = 1:size(bot_fills,1)
                for i_top = 1:size(top_fills,2)
                    if ~enforce_min_feat_size_func( periods(i_bot,i_top),...
                                                    top_fills(i_bot,i_top), bot_fills(i_bot,i_top), offsets(i_bot,i_top) )
                        directivities(i_bot,i_top) = 0;
                    end
                end
            end

            % for each top fill, pick bottom fill with highest
            % directivity
            [ chosen_dirs_vs_top, indx_highest_dir_per_topfill ] = max( directivities, [], 1 );
            % gotta use linear indexing
            indxs_vs_top               = sub2ind( size(directivities), indx_highest_dir_per_topfill, 1:size(directivities,2)  );

            % for each bottom fill, pick top fill
            % with highest directivity
            % directivity
            [ chosen_dirs_vs_bot, indx_highest_dir_per_botfill ] = max( directivities, [], 2 );
            % gotta use linear indexing
            indxs_vs_bot                   = sub2ind( size(directivities), 1:size(directivities,1), indx_highest_dir_per_botfill.' );

            switch fill_vs_top_bot_both

                case 'top'

                    % use top fill
                    chosen_angles       = angles( indxs_vs_top );
                    chosen_periods      = periods( indxs_vs_top );
                    chosen_offsets      = offsets( indxs_vs_top );
                    chosen_scatter_str  = scatter_strs( indxs_vs_top );
                    chosen_ks           = ks( indxs_vs_top );
                    chosen_top_fills    = top_fills( indxs_vs_top );
                    chosen_bot_fills    = bot_fills( indxs_vs_top );
                    chosen_dirs         = chosen_dirs_vs_top;

                case 'bottom'

                    % use bot fill
                    chosen_angles       = angles( indxs_vs_bot );
                    chosen_periods      = periods( indxs_vs_bot );
                    chosen_offsets      = offsets( indxs_vs_bot);
                    chosen_scatter_str  = scatter_strs( indxs_vs_bot );
                    chosen_ks           = ks( indxs_vs_bot );
                    chosen_top_fills    = top_fills( indxs_vs_bot );
                    chosen_bot_fills    = bot_fills( indxs_vs_bot );
                    chosen_dirs         = chosen_dirs_vs_bot;

                case 'both'

                    % append bot fill
                    chosen_angles       = [ angles( indxs_vs_top ), angles( indxs_vs_bot ) ];
                    chosen_periods      = [ periods( indxs_vs_top ), periods( indxs_vs_bot ) ];
                    chosen_offsets      = [ offsets( indxs_vs_top ), offsets( indxs_vs_bot) ];
                    chosen_scatter_str  = [ scatter_strs( indxs_vs_top ), scatter_strs( indxs_vs_bot ) ];
                    chosen_ks           = [ ks( indxs_vs_top ), ks( indxs_vs_bot ) ];
                    chosen_top_fills    = [ top_fills( indxs_vs_top ), top_fills( indxs_vs_bot ) ];
                    chosen_bot_fills    = [ bot_fills( indxs_vs_top ), bot_fills( indxs_vs_bot ) ];
                    chosen_dirs = [ chosen_dirs_vs_top, chosen_dirs_vs_bot.' ];

            end

            % sorting on scattering strength (might not be needed)
            [ chosen_scatter_str, indx_sort ] = sort( chosen_scatter_str );
            chosen_angles       = chosen_angles( indx_sort );
            chosen_periods      = chosen_periods( indx_sort );
            chosen_offsets      = chosen_offsets( indx_sort );
            chosen_ks           = chosen_ks( indx_sort );
            chosen_top_fills    = chosen_top_fills( indx_sort );
            chosen_bot_fills    = chosen_bot_fills( indx_sort );
            chosen_dirs         = chosen_dirs( indx_sort );
            
            % remove 0 directionality cells
            chosen_scatter_str  = chosen_scatter_str( chosen_dirs > 0 );
            chosen_angles       = chosen_angles( chosen_dirs > 0 );
            chosen_periods      = chosen_periods( chosen_dirs > 0 );
            chosen_offsets      = chosen_offsets( chosen_dirs > 0 );
            chosen_ks           = chosen_ks( chosen_dirs > 0 );
            chosen_top_fills    = chosen_top_fills( chosen_dirs > 0 );
            chosen_bot_fills    = chosen_bot_fills( chosen_dirs > 0 );
            chosen_dirs         = chosen_dirs( chosen_dirs > 0 );
            
            % interpolating here
            if interp_on
                n_interp            = 200;
                chosen_scatter_str  = interp1( 1:length(chosen_scatter_str), chosen_scatter_str, linspace(1, length(chosen_scatter_str), n_interp) );
                chosen_angles       = interp1( 1:length(chosen_angles), chosen_angles, linspace(1, length(chosen_angles), n_interp) );
                chosen_periods      = interp1( 1:length(chosen_periods), chosen_periods, linspace(1, length(chosen_periods), n_interp) );
                chosen_offsets      = interp1( 1:length(chosen_offsets), chosen_offsets, linspace(1, length(chosen_offsets), n_interp) );
                chosen_ks           = interp1( 1:length(chosen_ks), chosen_ks, linspace(1, length(chosen_ks), n_interp) );
                chosen_top_fills    = interp1( 1:length(chosen_top_fills), chosen_top_fills, linspace(1, length(chosen_top_fills), n_interp) );
                chosen_bot_fills    = interp1( 1:length(chosen_bot_fills), chosen_bot_fills, linspace(1, length(chosen_bot_fills), n_interp) );
                chosen_dirs         = interp1( 1:length(chosen_dirs), chosen_dirs, linspace(1, length(chosen_dirs), n_interp) );
            end
                
            % enforce min feature size
            obj.synthesized_design.use_min_feat_size = false;                   % default to false
            if exist( 'enforce_min_feat_size_func', 'var' )
                % loop through each cell and discard any that violate
                % feature size rules
                obj.synthesized_design.use_min_feat_size = true;

                indices_to_keep = [];
                for ii = 1:length( chosen_periods )
                    if enforce_min_feat_size_func( chosen_periods(ii), chosen_top_fills(ii), chosen_bot_fills(ii), chosen_offsets(ii) ) == true
                        indices_to_keep(end+1) = ii;
                    end
                end

                chosen_angles       = chosen_angles(indices_to_keep);
                chosen_periods      = chosen_periods(indices_to_keep);
                chosen_offsets      = chosen_offsets(indices_to_keep);
                chosen_scatter_str  = chosen_scatter_str(indices_to_keep);
                chosen_ks           = chosen_ks(indices_to_keep);
                chosen_top_fills    = chosen_top_fills(indices_to_keep);
                chosen_bot_fills    = chosen_bot_fills(indices_to_keep);
                chosen_dirs         = chosen_dirs(indices_to_keep);

            end
           
            if plot_on
                OPTS = struct('fig_size', [2400, 1200]);
                figure('Name', 'chosen_datapoints');
                % chosen angles
                subplot(3,3,1);
                plot( 1:length(chosen_angles), chosen_angles, '-o' );
                title('chosen datapoints, angles');
                xlim( [1, length(chosen_angles) ] );
                makeFigureNice( OPTS );
                % chosen periods
                subplot(3,3,2);
                plot( 1:length(chosen_periods), chosen_periods, '-o' );
                title('chosen datapoints, periods');
                makeFigureNice( OPTS );
                % chosen offsets
                subplot(3,3,3);
                plot( 1:length(chosen_offsets), chosen_offsets, '-o' );
                title('chosen datapoints, offsets');
                makeFigureNice( OPTS );
                % chosen scatter strength
                subplot(3,3,4);
                plot( 1:length(chosen_scatter_str), chosen_scatter_str, '-o' );
                title('chosen datapoints, scattering strengths');
                makeFigureNice( OPTS );
                % chosen k real
                subplot(3,3,5);
                plot( 1:length(chosen_ks), real(chosen_ks), '-o' );
                title('chosen datapoints, k real');
                makeFigureNice( OPTS );
                % chosen k imag
                subplot(3,3,6);
                plot( 1:length(chosen_ks), imag(chosen_ks), '-o' );
                title('chosen datapoints, k imaginary');
                makeFigureNice( OPTS );
                % chosen top fill
                subplot(3,3,7);
                plot( 1:length(chosen_top_fills), chosen_top_fills, '-o' );
                title('chosen datapoints, top fill');
                makeFigureNice( OPTS );
                % chosen bottom fill
                subplot(3,3,8);
                plot( 1:length(chosen_bot_fills), chosen_bot_fills, '-o' );
                title('chosen datapoints, bottom fill');
                makeFigureNice( OPTS );
                % chosen directivity
                subplot(3,3,9);
                plot( 1:length(chosen_dirs), 10*log10(chosen_dirs), '-o' );
                title('chosen datapoints, directivity (dB)');
                makeFigureNice( OPTS );
            end

            % save the chosen cells
            obj.chosen_cells = struct( 'chosen_angles', chosen_angles, ...
                                    'chosen_periods', chosen_periods, ...
                                    'chosen_offsets', chosen_offsets, ...
                                    'chosen_scatter_str', chosen_scatter_str, ...
                                    'chosen_ks', chosen_ks, ...
                                    'chosen_top_fills', chosen_top_fills, ...
                                    'chosen_bot_fills', chosen_bot_fills, ...
                                    'chosen_dirs', chosen_dirs );
                                    
        end % end choose_unit_cells()
        
       
        
        function [ obj, best_alpha_power ] = optimize_start_alpha( obj, xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, desired_field )
            % Optimizes the starting alpha, based on predicted field overlap
            %
            % Inputs:
            %     desired_field
            %         type: double, array
            %         desc: desired field vs. xvec
            
%             % pick the alphas to try:
%             alpha_power_max = log10( min( scatter_strs_high_dir ) );
%             alpha_powers    = linspace( -6, alpha_power_max, 100 );
%             alpha_try       = 10.^( alpha_powers );
            
%             % save max overlap, dimensions overlap vs. alpha
%             max_overlap = zeros( size( alpha_try ) );
            
%             % brute force optimize
%             for i_alpha = 1:length(alpha_powers)
%                 
%                 max_overlap(i_alpha) = obj.predict_overlap_for_optimization( ...
%                                               xvec, alpha_des, high_dirs, ...
%                                               bot_fills_high_dir, top_fills_high_dir, ...
%                                               offsets_high_dir, periods_high_dir, ...
%                                               angles_high_dir, scatter_strs_high_dir, ...
%                                               k_high_dir, desired_field, alpha_powers(i_alpha) );
%                 
%             end     % end for i_alpha = 1:length(alpha_try)
           
            % using matlab's fminsearch
            f = @(alpha_power)( 1 - obj.predict_overlap_for_optimization( ...
                                              xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, desired_field, alpha_power ) );
            
            % options
%             opts = optimset( 'Display', 'iter' );
                                          
            % for debuggging
%             [ best_alpha_power, best_val] = fminsearch( f, log10(min( scatter_strs_high_dir )), opts );
            
            best_alpha_power = fminsearch( f, log10(min( chosen_scatter_strs )) );
            
%             best_overlap_fminsearch = 1 - best_val;
            
%             % DEBUG plot max overlap vs. alpha
%             figure;
%             plot( alpha_powers, max_overlap, '-o' );
%             xlabel('alpha power'); ylabel('max overlap');
%             title('DEBUG max overlap vs. start alpha (in power)');
%             makeFigureNice();
            
            % now save the final design with the best predicted overlap
            obj = obj.pick_final_datapoints(  xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, 10.^(best_alpha_power) );
                                          
        end     % end function optimize_start_alpha()
        

        function [ max_overlap, field_shape_prediction ] = predict_overlap_for_optimization( ...
                                              obj, xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, desired_field, start_alpha_power )
            % merit function used in optimize_start_alpha for predicting
            % overlap
            %
            % Inputs:
            %   start_alpha_power:
            %       type: double, scalar
            %       desc: alpha to try = 10^(start_alpha_power)
            %
            % Outputs:
            %     max_overlap
            %     field_shape_prediction
            %         type: double, array
            %         desc: predicted field shape vs. x, mostly for debugging
                      
            % pick the datapoints
            [ obj, synthesized_design ] = obj.pick_final_datapoints(  xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, 10.^(start_alpha_power) );
                                          
            % convert scatter str vs. cell into scatter str vs. x
            xvec                = xvec - xvec(1);
            scatter_str_vs_x    = zeros( size( xvec ) );
            end_pos_vs_cell     = cumsum( synthesized_design.period );
            cur_cell            = 1;
            for ii = 1:length(xvec)

                scatter_str_vs_x(ii) = synthesized_design.scatter_str(cur_cell);

                if xvec(ii) > end_pos_vs_cell( cur_cell ) && cur_cell < length(end_pos_vs_cell)
                    % move to next cell
                    cur_cell = cur_cell + 1;
                end

            end

            % integrate the picked scatter strengths to get the predicted
            % field shape
            
            [obj, field_shape_prediction] = obj.predict_field_shape( xvec, scatter_str_vs_x );
            
            % calculate overlap with desired field
            [ ~, overlap ]  = obj.xcorr_normalized( desired_field, field_shape_prediction );
            max_overlap     = max( abs(overlap) );
            
        end     % end function predict_overlap_for_optimization()
  
        
        function [ obj, synthesized_design ] = pick_final_datapoints( obj, xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, start_alpha_des )
            % Picks the final datapoints (cells) that make up the grating
            %
            % Inputs:
            %     xvec
            %         type: double, array
            %         desc: coordinates of desired field profile
            %     alpha_des
            %         type: double, array
            %         desc: desired scattering strength vs. x
            %     high_dirs
            %         type: double, array
            %         desc: chosen highest directivity points
            %     bot_fills_high_dir
            %         type: double, array
            %         desc: bottom fill value of chosen highest dir. points
            %     top_fills_high_dir
            %         type: double, array
            %         desc: top fill value of chosen highest dir. points
            %     offsets_high_dir
            %         type: double, array
            %         desc: offset values of chosen highest dir. points
            %     periods_high_dir
            %         type: double, array
            %         desc: periods of chosen highest dir. points
            %     angles_high_dir
            %         type: double, array
            %         desc: angle of chosen highest dir. points
            %     scatter_strs_high_dir
            %         type: double, array
            %         desc: scattering strength (alpha) of chosen highest
            %               dir. points
            %     k_high_dir
            %         type: double, array
            %         desc: k of chosen highest dir. points
            %     start_alpha_des
            %         type: double, scalar
            %         desc: OPTIONAL, user defined starting alpha
            %
            % Outputs:
            %   synthesized_design:
            %       type: struct
            %       desc: a struct that holds the synthesized cells
            
            % default to weakest cell
            if ~exist( 'start_alpha_des', 'var' )
                start_alpha_des = min(chosen_scatter_strs); %1e-5;
            end
            
            % now match these data points to the desired alpha
            % starting point
            [~, indx_max_alpha] = max( alpha_des );
            [~, indx_x]         = min(abs( alpha_des(1:indx_max_alpha) - start_alpha_des ) );
            cur_x               = xvec(indx_x);
            
            % final synthesized variables
            synthesized_design.dir                  = [];
            synthesized_design.bot_fill             = [];
            synthesized_design.top_bot_fill_ratio   = [];
            synthesized_design.top_fill             = [];
            synthesized_design.period               = [];
            synthesized_design.offset               = [];
            synthesized_design.angles               = [];
            synthesized_design.scatter_str          = [];
            synthesized_design.k                    = [];
            synthesized_design.GC                   = {};
            synthesized_design.des_scatter          = [];
            
            
            % flag for switching to using max scattering strength
            saturate_scatter_str_to_max = false;
 
            ii = 1;
            while cur_x < xvec(end)
                % build grating one cell at a time
                
                % pick design with scattering strength closest to desired
                % alpha
                des_scatter = alpha_des(indx_x);                            % desired alpha
                if ii > 1
                    if des_scatter  > max( chosen_scatter_strs ) || des_scatter < synthesized_design.des_scatter(ii-1)
                        % desired scattering strength too high, gotta saturate
                        saturate_scatter_str_to_max = true;
                    end
                end
                if ~saturate_scatter_str_to_max
                    [~, indx_closest_scatter]   = min( abs(chosen_scatter_strs - des_scatter) );          % index of closest scatter design 
                else
                    [~, indx_closest_scatter]   = max( chosen_scatter_strs );                             % saturate to max
                end
                
                % save parameters
                synthesized_design.dir(ii)                   = chosen_dirs( indx_closest_scatter );
                synthesized_design.bot_fill(ii)              = chosen_bot_fills( indx_closest_scatter );
                synthesized_design.top_fill(ii)              = chosen_top_fills( indx_closest_scatter );
                synthesized_design.top_bot_fill_ratio(ii)    = chosen_top_fills( indx_closest_scatter ) ./ chosen_bot_fills( indx_closest_scatter );
                synthesized_design.offset(ii)                = chosen_offsets( indx_closest_scatter );
                synthesized_design.period(ii)                = chosen_periods( indx_closest_scatter );
                synthesized_design.angles(ii)                = chosen_angles( indx_closest_scatter );
                synthesized_design.scatter_str(ii)           = chosen_scatter_strs( indx_closest_scatter );
                synthesized_design.k(ii)                     = chosen_ks( indx_closest_scatter );
                synthesized_design.des_scatter(ii)           = des_scatter;
                
                synthesized_design.GC{ii} = obj.h_makeGratingCell(    ...
                               obj.discretization, ...
                               obj.background_index, ...
                               obj.y_domain_size, ...
                               synthesized_design.period(ii), ...
                               synthesized_design.top_fill(ii), ...
                               synthesized_design.bot_fill(ii), ...
                               synthesized_design.offset(ii)/synthesized_design.period(ii) );
                           
                % move onto next
                cur_x       = cur_x + synthesized_design.period(ii);
                [~, indx_x] = min( abs(xvec - cur_x) );
                cur_x       = xvec( indx_x );
                ii          = ii + 1;
                
            end     % end for ii = 1:ncells
            
        end     % end pick_final_datapoints()
              
        
        function [ obj, top_fills, bot_fills, directivities, ...
                   angles, periods, offsets, scatter_strs, ks, rad_power_ratio ] ...
                   = pick_design_quadrant( obj, input_wg_type, fill_top_override, fill_bot_override )
            % Picks out the design quadrant of the top/bot fill space that
            % corresponds to the desired input waveguide type
            %
            % Inputs:
            %   input_wg_type
            %       type: string
            %       desc: 'top' 'bottom' 'full' 'none'
            %   fill_top_override
            %       type: scalar
            %       desc: OPTIONAL - overrides the center of mass
            %   fill_bot_override
            %       type: scalar
            %       desc: OPTIONAL - overrides the center of mass 
            
            % meshgrid the fills (dimensions are bot vs. top, I believe)
            [ top_fills_mesh, bot_fills_mesh ] = meshgrid( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots );

            % find where the max scattering strength is
            if nargin < 3
                % new ver. where i find center of mass
                r       = [ top_fills_mesh(:), bot_fills_mesh(:) ];   % dimensions [ top fill vs. mass, bot fill vs. mass ]
                mass    = obj.sweep_variables.scatter_str_vs_fills(:); 
                R       = sum( r .* [ mass, mass ], 1 )./sum(mass(:));
                [~, i_max_top_scat] = min(abs( obj.sweep_variables.fill_tops - R(1) ) );
                [~, i_max_bot_scat] = min(abs( obj.sweep_variables.fill_bots - R(2) ) );
            else
                % override center of mass
                [~, i_max_top_scat] = min(abs( obj.sweep_variables.fill_tops - fill_top_override ) );
                [~, i_max_bot_scat] = min(abs( obj.sweep_variables.fill_bots - fill_bot_override ) );
            end
            
            % Select subset of domain to use
            if strcmp( input_wg_type, 'bottom' ) == true
                % input waveguide is the bottom layer
            
                % narrow down design space (take top left quadrant)
                % recall fills are sorted in descending order
                top_fills       = top_fills_mesh( 1:i_max_bot_scat, i_max_top_scat:end );
                bot_fills       = bot_fills_mesh( 1:i_max_bot_scat, i_max_top_scat:end );
                directivities   = obj.sweep_variables.directivities_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                angles          = obj.sweep_variables.angles_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );      
                periods         = obj.sweep_variables.periods_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                offsets         = obj.sweep_variables.offsets_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                ks              = obj.sweep_variables.k_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                try
                    rad_power_ratio = obj.sweep_variables.prad_pin_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                catch
                end
                
            elseif strcmp( input_wg_type, 'full' ) == true
                % input waveguide is a full waveguide
                
                % narrow down design space (take top right quadrant)
                top_fills       = top_fills_mesh( 1:i_max_bot_scat, 1:i_max_top_scat );
                bot_fills       = bot_fills_mesh( 1:i_max_bot_scat, 1:i_max_top_scat );
                directivities   = obj.sweep_variables.directivities_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                angles          = obj.sweep_variables.angles_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );      
                periods         = obj.sweep_variables.periods_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                offsets         = obj.sweep_variables.offsets_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                ks              = obj.sweep_variables.k_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                try
                    rad_power_ratio = obj.sweep_variables.prad_pin_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                catch
                end
                    
            elseif strcmp( input_wg_type, 'top' ) == true
                % input waveguide is a full waveguide
                
                % narrow down design space (take bottom right quadrant)
                top_fills       = top_fills_mesh( i_max_bot_scat:end, 1:i_max_top_scat );
                bot_fills       = bot_fills_mesh( i_max_bot_scat:end, 1:i_max_top_scat );
                directivities   = obj.sweep_variables.directivities_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                angles          = obj.sweep_variables.angles_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );      
                periods         = obj.sweep_variables.periods_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                offsets         = obj.sweep_variables.offsets_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                ks              = obj.sweep_variables.k_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                try
                    rad_power_ratio = obj.sweep_variables.prad_pin_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                catch
                end
                
            elseif strcmp( input_wg_type, 'none' ) == true
                % input waveguide is cladding
                
                % narrow down design space (take bottom left quadrant)
                top_fills       = top_fills_mesh( i_max_bot_scat:end, i_max_top_scat:end );
                bot_fills       = bot_fills_mesh( i_max_bot_scat:end, i_max_top_scat:end );
                directivities   = obj.sweep_variables.directivities_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                angles          = obj.sweep_variables.angles_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );      
                periods         = obj.sweep_variables.periods_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                offsets         = obj.sweep_variables.offsets_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                ks              = obj.sweep_variables.k_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                try
                    rad_power_ratio = obj.sweep_variables.prad_pin_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                catch
                end
                
            end     % end if strcmp( input_wg_type, 'bottom' )
            
        end     % end pick_design_quadrant()

        function [] = plot_chosen_cells( obj )
            % plotting chosen cells from stage 2, prior to concatenation
            
            N_cells = length(obj.chosen_cells.chosen_angles);
            OPTS = struct('fig_size', [2400, 1200]);
            figure('Name', 'chosen_datapoints');
            % chosen angles
            subplot(3,3,1);
            plot( 1:N_cells, obj.chosen_cells.chosen_angles, '-o' );
            title('chosen datapoints, angles');
            xlim( [1, length(chosen_angles) ] );
            makeFigureNice( OPTS );
            % chosen periods
            subplot(3,3,2);
            plot( 1:N_cells, obj.chosen_cells.chosen_periods, '-o' );
            title('chosen datapoints, periods');
            makeFigureNice( OPTS );
            % chosen offsets
            subplot(3,3,3);
            plot( 1:N_cells, obj.chosen_cells.chosen_offsets, '-o' );
            title('chosen datapoints, offsets');
            makeFigureNice( OPTS );
            % chosen scatter strength
            subplot(3,3,4);
            plot( 1:N_cells, obj.chosen_cells.chosen_scatter_str, '-o' );
            title('chosen datapoints, scattering strengths');
            makeFigureNice( OPTS );
            % chosen k real
            subplot(3,3,5);
            plot( 1:N_cells, real(obj.chosen_cells.chosen_ks), '-o' );
            title('chosen datapoints, k real');
            makeFigureNice( OPTS );
            % chosen k imag
            subplot(3,3,6);
            plot( 1:N_cells, imag(obj.chosen_cells.chosen_ks), '-o' );
            title('chosen datapoints, k imaginary');
            makeFigureNice( OPTS );
            % chosen top fill
            subplot(3,3,7);
            plot( 1:N_cells, obj.chosen_cells.chosen_top_fills, '-o' );
            title('chosen datapoints, top fill');
            makeFigureNice( OPTS );
            % chosen bottom fill
            subplot(3,3,8);
            plot( 1:N_cells, obj.chosen_cells.chosen_bot_fills, '-o' );
            title('chosen datapoints, bottom fill');
            makeFigureNice( OPTS );
            % chosen directivity
            subplot(3,3,9);
            plot( 1:N_cells, 10*log10(obj.chosen_cells.chosen_dirs), '-o' );
            title('chosen datapoints, directivity (dB)');
            makeFigureNice( OPTS );
            
        end
        
    end     % End methods section
    
end     % end class definition
















































