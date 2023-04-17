classdef c_synthGrating
% Synthesizes a grating that outputs a desired field profile
%
% Authors: bohan zhang
%
% Prerequisites/dependencies
%   - c_gratingCell, or some child class of c_gratingCell
%   - the utility folder ?
%
%
%   The user should define their own custom grating unit cell
%   drawing function.
%   HOWEVER, this function MUST have the following inputs and outputs, IN
%   ORDER:
%       function GC = your_makeGratingCell_function( dxy, background_index, y_domain_size, period, fill )
%           % makes and returns a c_gratingCell object (or subclass of one)
% 
%       outputs:
%           GC
%               type: c_gratingCell object
%               desc: single level grating cell object
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
%
%
% Examples:
%
% Example code that makes a synthesis object and generates the design space
%
%   % dependencies
%   addpath( genpath( 'C:\Users\bz\git\grating_synthesis' ) );                  % path on plab desktop to grating synth codes
% 
%   % inputs
%   discretization      = 10;
%   units               = 'nm';
%   lambda              = 1550;
%   background_index    = 1.0;
%   y_domain_size       = 5000;
%   optimal_angle       = -15;
%   data_notes          = 'whatever';
% 
%   % grating cell making function handle
%   h_makeGratingCell = @f_simple_onelevel_cell;
% 
%   % make synthesis object
%   synth_obj = c_synthGrating( 'discretization',    discretization, ...
%                               'units',             units,   ...
%                               'lambda',            lambda, ...
%                               'background_index',  background_index,    ...
%                               'y_domain_size',     y_domain_size, ...
%                               'optimal_angle',     optimal_angle, ...
%                               'data_notes',        data_notes, ...
%                               'h_makeGratingCell', h_makeGratingCell ...
%                               );
%                         
%   % generate design space
%   synth_obj   = synth_obj.generate_design_space();
%   sweep_vars  = synth_obj.sweep_variables;



    properties

        discretization;     % dx and dy
        units;              % units, verbose, 'm' or 'mm', or 'um', or 'nm'
                            % has fields 'name' and 'scale'
        lambda;             % center wavelength
        k0;                 % center wavenumber
        background_index;   % background index
        y_domain_size;      % transverse domain size (y coordinate aka vertical dimension)
        inputs;             % saves input settings for user reference
        start_time;         % time when object was created, 'YEAR-month-day hour-min-sec'
        data_notes;         % verbose notes of what current sweep is doing
                            
        h_makeGratingCell;  % handle to the grating cell making function
        
        % parameters to optimize for
        coupling_direction;     % either 'up' or 'down
        optimal_angle;          % angle to optimize for, deviation from the normal, in deg.
        
        % struct that holds resulting variables from sweep
        % has the following fields:
        % fill_ratios_to_sweep
        % directivities_vs_fill 
        % angles_vs_fill
        % scatter_str_vs_fill
        % periods_vs_fill
        % k_vs_fill
        % GC_vs_fill
        % overlap_predict_vs_fill - if uniform design is run
        sweep_variables;
        
        % struct that stores final design parameters
        % currently stores the following:
        % GC_synth
        % scatter_str_synth
        % fill_synth
        % N
        % 
        synthesized_design;
         
    end

    
    
    methods
        
        function obj = c_synthGrating(varargin)
            % Constructor
            % See top comments for input documentation
              
            % inputs and defaults
            inputs = {  'discretization',   'none', ...
                        'lambda',           'none', ...
                        'background_index', 1.0,    ...
                        'y_domain_size',    'none', ...
                        'optimal_angle',    'none', ...
                        'data_notes',       '', ...
                        'coupling_direction',   'none', ...
                        'h_makeGratingCell', @makeGratingCell ...
                     }; 
            obj.inputs = inputs;
            
            % parse inputs
            p = f_parse_varargin( inputs, varargin{:} );

            % save starting time
            obj.start_time = datestr( datetime('now'), 'yyyy_mm_dd_HH_MM_SS_' );

            % set other properties
            obj.discretization      = p.discretization;
            obj.lambda              = p.lambda;
            obj.k0                  = 2*pi/p.lambda;
            obj.background_index    = p.background_index;
            obj.y_domain_size       = p.y_domain_size;
            obj.optimal_angle       = p.optimal_angle;

            if strcmp( p.coupling_direction, 'up') || strcmp( p.coupling_direction, 'down') 
                % set coupling direction
                obj.coupling_direction = p.coupling_direction;
            else
                error('Error: input ''coupling_direction'' is not valid. Valid entries are ''up'' or ''down''. You entered ''%s''', p.coupling_direction);
            end

            % user notes/documentation
            obj.data_notes      = p.data_notes;
                
            % set handle to grating cell making function
            obj.h_makeGratingCell = p.h_makeGratingCell;

        end     % end constructor()
               
        
        function save_to_struct(obj, filename)
            % Saves all current properties of this object to a structure,
            % and then to a .mat file
            sweep_obj = obj.convertObjToStruct();
            save(filename, 'sweep_obj');
        end
        
        
        function obj_as_struct = convert_obj_to_struct(obj)
            % converts the current object to a struct that holds the
            % object's properties
            
            props = properties(obj);
     
            obj_as_struct = struct();
            for p = 1:numel(props)
                if strcmp( props{p}, 'h_makeGratingCell' )
                    % convert function handle to string
                    obj_as_struct.(props{p}) = func2str( obj.(props{p}) );
                else
                    obj_as_struct.(props{p}) = obj.(props{p});
                end
            end
            
        end
            
        function [obj, u] = fiber_mode_gaussian(obj, w0, zvec, xvec, theta, d0, nclad)
            % somewhat adapted from Cale's code
            %
            % Generate Gaussian-beam mode profile at a plane through y = d0 at angle of
            % theta (in degrees)
            %
            % (using H.A. Haus, Waves & Fields, Chapter 5)
            % CMG November 21st, 2014
            %
            %
            % Inputs:   
            %   w0  
            %       type: double, scalar
            %       desc: 1/e beam RADIUS at waist, in units 'units'
            %           
            %   zvec
            %       type: double, array
            %       desc: out of plane coordinates... currently not used
            %
            %   xvec
            %       type: double, array
            %       desc: coordinates along direction of propagation of
            %             grating, in units 'units'
            %           
            %   theta 
            %       type: double, scalar
            %       desc: angle from normal in degrees
            %
            %   d0  
            %       type: double, scalar
            %       desc: distance from beam waist to slice
            %
            %   nclad 
            %       type: double, scalar
            %       desc: cladding index
            %
            %
            % Outputs: 
            %   u
            %       type: double, array
            %       desc: returned slice of gaussian beam, normalized to total
            %             power


            % Constants, in units of meters
            lambda  = obj.lambda/nclad; % * obj.units.scale / nclad;     % wavelength in cladding, units m
            k0      = 2*pi/lambda;                              % 1/m
            
            % Convert to radians
            theta = (pi/180)*theta;
            
            % Scale coordinates% units m
            yvec = xvec;
            
            % dimensions in frame of gaussian beam optical axis 
            xprime = xvec.*cos(-theta) + d0*sin(-theta);
            zprime = -xvec.*sin(-theta) + d0*cos(-theta);

            % b (confocal parameters) is used instead of z0 so that z0 = -1j.*b removes the singularity of the solution on the real z axis (see Haus pg 109)
            b = k0*w0^2/2;                                                                                   

            % Equation (5.2) in Haus
            u00_slice =   1j .* sqrt(k0*b/pi) .* ( 1./(zprime + 1j.*b) ).*...
                exp( -1j.*k0.*( xprime.^2 )./( 2*(zprime + 1j.*b) ) );     
            
            % normalize the slice to intensity
            dx          = obj.discretization; % * obj.units.scale;                 % disc. in m
            u00_slice   = u00_slice/sqrt( dx * sum( abs( u00_slice ).^2 ) );
            
            % return and save data
            u       = u00_slice;

        end     % end fiber_mode_gaussian()
        
        
        function obj = generate_design_space( obj, fill_ratios_to_sweep )
            % sweep fills, optimize period for a single output angle
            %
            % Inputs:
            %   fill_ratios_to_sweep
            %       type: double, array
            %       desc: OPTIONAL fill ratios to sweep
            
            % default fills
            if nargin < 2
                fill_ratios_to_sweep = fliplr( 0.02:0.02:0.98 );
            end
            
            % sort fills so they are in descending order
            fill_ratios_to_sweep = sort( fill_ratios_to_sweep, 'descend' );
            
            % make waveguide cell
            waveguide = obj.h_makeGratingCell( obj.discretization, .....
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               2*obj.discretization, ...
                                               1.0 );
            
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
            waveguide_k = waveguide.k;                                      % units of rad/'units'                    
            
            % calculate analytical period which would approximately phase
            % match to desired output angle
            k0              = obj.background_index * obj.k0;
            kx              = k0 * sin( (pi/180) * obj.optimal_angle );
            guess_period    = 2*pi/(waveguide_k - kx);                              % units of 'units'
            
            % snap period to discretization
            guess_period    = obj.discretization * round(guess_period/obj.discretization);

            % ugh this is really annoying but i have to extend the
            % waveguide's e z overlap
            [ waveguide, e_z_overlap_ext ]  = waveguide.stitch_E_field( waveguide.Phi, real(waveguide.k), round(guess_period/waveguide.domain_size(2)) );
            waveguide.E_z_for_overlap       = e_z_overlap_ext;
            
            % initially start with waveguide GC
            guess_GC = waveguide;
            
            % set grating solver settings
            num_modes   = 5;
            BC          = 0;                                                % 0 = PEC
            pml_options = [1, 100, 20, 2]; 
            OPTS        = struct( 'mode_to_overlap', e_z_overlap_ext );
     
            % initialize saving variables
            obj.sweep_variables.fill_ratios_to_sweep    = fill_ratios_to_sweep;
            obj.sweep_variables.directivities_vs_fill   = zeros( size(fill_ratios_to_sweep) );
            obj.sweep_variables.angles_vs_fill          = zeros( size(fill_ratios_to_sweep) );
            obj.sweep_variables.scatter_str_vs_fill     = zeros( size(fill_ratios_to_sweep) );
            obj.sweep_variables.periods_vs_fill         = zeros( size(fill_ratios_to_sweep) );
            obj.sweep_variables.k_vs_fill               = zeros( size(fill_ratios_to_sweep) );
            obj.sweep_variables.GC_vs_fill              = cell( size(fill_ratios_to_sweep) );
            
            % for each fill, optimize the period to achieve closest angle
            tic;
            for ii = 1:length(fill_ratios_to_sweep)
               
                fprintf('Sweeping fill ratio %i of %i\n', ii, length(fill_ratios_to_sweep));
                
                % simulate the grating, get the angle
                GC = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               guess_period, ...
                                               fill_ratios_to_sweep(ii) );
                GC = GC.runSimulation( num_modes, BC, pml_options, obj.k0, guessk, OPTS );
                
                % init saving variables
                angles_vs_period    = [];
                k_vs_period         = [];
                GC_vs_period        = {};
                periods             = [];
                
                % update for next iteration
                periods(1)              = GC.domain_size(2);
                GC_vs_period{1}         = GC;
                k_vs_period(1)          = GC.k;
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    angles_vs_period( 1 ) = GC.max_angle_up;
                else
                    % coupling direction is downwards
                    angles_vs_period( 1 ) = GC.max_angle_down;
                end
                guessk                  = GC.k;
                OPTS.mode_to_overlap    = GC.E_z_for_overlap;
                
                % decide whether to sweep larger or smaller periods
                % based on the angle
                if angles_vs_period(1) > obj.optimal_angle
                    % only sweep smaller periods
                    delta_period = -obj.discretization;
                else
                    % only sweep larger periods
                    delta_period = obj.discretization;
                end
            
                i_period = 2;
                while true
                   
                    fprintf('Period iteration %i\n', i_period );
                    
                    % update period
                    guess_period = guess_period + delta_period;

                    % make grating cell
                    GC = obj.h_makeGratingCell( obj.discretization, ...
                                                   obj.background_index, ...
                                                   obj.y_domain_size, ...
                                                   guess_period, ...
                                                   fill_ratios_to_sweep(ii) );

                    % run sim
                    GC = GC.runSimulation( num_modes, BC, pml_options, obj.k0, guessk, OPTS );

                    if strcmp( obj.coupling_direction, 'up' )
                        % coupling direction is upwards
                        angles_vs_period( i_period ) = GC.max_angle_up;
                    else
                        % coupling direction is downwards
                        angles_vs_period( i_period ) = GC.max_angle_down;
                    end

                    % update for next iteration
                    periods(i_period)           = guess_period;
                    GC_vs_period{i_period}      = GC;
                    k_vs_period(i_period)       = GC.k;
                    guessk                      = GC.k;
                    OPTS.mode_to_overlap        = GC.E_z_for_overlap;
                    
                    % check for exit condition (if error in angle gets worse)
                    cur_angle_err   = abs( angles_vs_period( i_period ) - obj.optimal_angle );
                    prev_angle_err  = abs( angles_vs_period( i_period-1 ) - obj.optimal_angle );
                    if cur_angle_err > prev_angle_err
                        % optimization over, break
                        break;
                    end
                
                    i_period = i_period + 1;
                    
                    toc;
                    
                end     % end period sweep

                % pick best period
                [angle_error, indx_best_period] = min( abs( obj.optimal_angle - angles_vs_period ) );
                best_GC                         = GC_vs_period{ indx_best_period };
                
                % coupling direction is assumed to be upwards
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    obj.sweep_variables.directivities_vs_fill( ii )    = best_GC.directivity;
                    obj.sweep_variables.angles_vs_fill( ii )           = best_GC.max_angle_up;
                    obj.sweep_variables.scatter_str_vs_fill( ii )      = best_GC.alpha_up;
%                     obj.sweep_variables.scatter_str_vs_fill_srad( ii )      = best_GC.alpha_up_from_srad;
                else
                    % coupling direction is downwards
                    obj.sweep_variables.directivities_vs_fill( ii )    = 1./best_GC.directivity;
                    obj.sweep_variables.angles_vs_fill( ii )           = best_GC.max_angle_down;
                    obj.sweep_variables.scatter_str_vs_fill( ii )      = best_GC.alpha_down;
                end
                obj.sweep_variables.periods_vs_fill( ii )          = best_GC.domain_size(2);
                obj.sweep_variables.k_vs_fill( ii )                = best_GC.k;
                
                % going to save this without the field to save memory
                obj.sweep_variables.GC_vs_fill{ ii } = obj.h_makeGratingCell( obj.discretization, ...
                                                   obj.background_index, ...
                                                   best_GC.domain_size(1), ...
                                                   best_GC.domain_size(2), ...
                                                   fill_ratios_to_sweep(ii) );
                
                % update guess period for next iteration
                guess_period = best_GC.domain_size(2);
                
                toc;
                
            end     % end for ii = 1:length(fill_ratios_to_sweep)
            
        end     % end generateDesignSpace()
        
       
        
        function obj = generate_final_design_uniform( obj, desired_field, input_wg_type, enforce_min_feat_size_func )
            % Synthesizes a uniform grating that has amplitude matching field_profile
            %
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
            %       desc: OPTIONAL user defined function for enforcing minimum feature sizes
            %               currently this function must take 2 args - ( period, fill )  
            
            % first choose which half of the design space to use
            [ obj, chosen_fills, chosen_periods, ...
                    chosen_angles, chosen_scatter_str, ...
                    chosen_ks ] = obj.pick_design_half( input_wg_type );
            
            % enforce min feature size
            obj.synthesized_design.use_min_feat_size = false;                   % default to false
            if exist( 'enforce_min_feat_size_func', 'var' )
                obj.synthesized_design.use_min_feat_size = true;
                
                [ obj, chosen_fills, chosen_periods, ... 
                chosen_angles, chosen_scatter_str, chosen_ks ]...
                = enforce_min_feat_size( obj, enforce_min_feat_size_func, ... 
                                            chosen_fills, ... 
                                            chosen_periods, ... 
                                            chosen_angles, ... 
                                            chosen_scatter_str, ... 
                                            chosen_ks );
            end
                
            % init saving vars
            overlap_predict_vs_fill = zeros( size( chosen_fills ) );
            
            % for each fill
            for i_fill = 1:length(chosen_fills)
                
                % grab current alpha
                cur_alpha = chosen_scatter_str(i_fill);
                
                % calc length of grating to scatter 99% of light (1/e^4 power)
                grat_len = 2/cur_alpha;     % in 'units'
                
                % gen x coords
                xvec = 0 : obj.discretization : grat_len - obj.discretization;
                
                % calculate the predicted field shape, which is simply an exponential decay
                field_shape_prediction = exp( -cur_alpha .* xvec );
                
                % calculate overlap with desired field (simple xcorr)
                [ ~, overlap ]  = obj.xcorr_normalized( desired_field, field_shape_prediction );
                overlap_predict_vs_fill(i_fill) = max( abs(overlap) );
                
            end
            
            % pick best overlap and synthesize final design with it
            [~, indx_best_overlap]  = max(abs(overlap_predict_vs_fill));
            final_grat_len          = 2 ./ chosen_scatter_str( indx_best_overlap );
            xvec                    = 0 : obj.discretization : final_grat_len - obj.discretization;
            alpha_des               = chosen_scatter_str( indx_best_overlap ) .* ones( size( xvec ) );
            
            [ obj, synthesized_design ] = obj.pick_final_datapoints( ...
                                      xvec, alpha_des, ...
                                      chosen_fills( indx_best_overlap ), ...
                                      chosen_periods( indx_best_overlap ), ...
                                      chosen_angles( indx_best_overlap ), ...
                                      chosen_scatter_str( indx_best_overlap ), ...
                                      chosen_ks( indx_best_overlap ) );
                                  
            obj.synthesized_design = catstruct( obj.synthesized_design, synthesized_design );   % combine structs
            obj.synthesized_design.input_wg_type = input_wg_type;
            
            % build final index distribution
            obj = obj.build_final_index();
            
        end     % end generate_final_design_uniform()
        
                
        function obj = generate_final_design_uniform_gaussian( obj, MFD, input_wg_type, enforce_min_feat_size_func )
            % Synthesizes an apodized grating that radiates desired Gaussian field profile
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter of gaussian, defined as 1/e^2 width of intensity
            %   input_wg_type
            %       type: string
            %       desc: 'full', 'bottom' or 'none' are currently supported
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: user defined function for enforcing minimum feature sizes
            %               currently this function must take 2 args - ( period, fill )   
            %             i guess this is optional
            
            % generate a fiber gaussian mode
            [ obj, field_profile ] = obj.make_gaussian_profile( MFD );
            
            % generate final design, apodized
            if exist( 'enforce_min_feat_size_func', 'var' )
                obj = obj.generate_final_design_uniform( field_profile, input_wg_type, enforce_min_feat_size_func );
            else
                obj = obj.generate_final_design_uniform( field_profile, input_wg_type );
            end
            
            % save MFD
            obj.synthesized_design.MFD = MFD;
            
        end     % end generate_final_design_uniform_gaussian()
        
        
        function obj = generate_final_design_apodized( obj, desired_field, input_wg_type, enforce_min_feat_size_func )
            % Synthesizes an apodized grating that radiates desired field profile
            %
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
            %               currently this function must take 2 args - ( period, fill )
            
            % calculate desired alpha
            [ obj, xvec, alpha_des, desired_field ] = obj.calculate_desired_scattering( desired_field );
            
            % first choose which half of the design space to use
            [ obj, chosen_fills, chosen_periods, ...
                    chosen_angles, chosen_scatter_str, ...
                    chosen_ks ] = obj.pick_design_half( input_wg_type );
            
            % enforce min feature size
            obj.synthesized_design.use_min_feat_size = false;                   % default to false
            if nargin > 3 %exist( 'enforce_min_feat_size_func', 'var' )
                % loop through each cell and discard any that violate
                % feature size rules
                obj.synthesized_design.use_min_feat_size = true;
                
                [ obj, chosen_fills, chosen_periods, ... 
                chosen_angles, chosen_scatter_str, chosen_ks ]...
                = enforce_min_feat_size( obj, enforce_min_feat_size_func, ... 
                                            chosen_fills, ... 
                                            chosen_periods, ... 
                                            chosen_angles, ... 
                                            chosen_scatter_str, ... 
                                            chosen_ks );
                 
            end
            
            % optimize for best starting alpha
            [ obj, best_alpha_power ] = obj.optimize_start_alpha( xvec, alpha_des,  ...
                                                chosen_fills, ...
                                                chosen_periods, ...
                                                chosen_angles, ...
                                                chosen_scatter_str, ...
                                                chosen_ks, ...
                                                desired_field );
            
            % synthesize final design
            [ obj, synthesized_design ] = obj.pick_final_datapoints( xvec, alpha_des, ...
                                      chosen_fills, ...
                                      chosen_periods, ...
                                      chosen_angles, ...
                                      chosen_scatter_str, ...
                                      chosen_ks, ...
                                      10.^(best_alpha_power) );
            obj.synthesized_design = catstruct( obj.synthesized_design, synthesized_design );   % combine structs
            obj.synthesized_design.input_wg_type = input_wg_type;
            
            % build final index distribution
            obj = obj.build_final_index();
            
        end     % end generate_final_design_apodized()
        
        
        function obj = generate_final_design_apodized_gaussian( obj, MFD, input_wg_type, enforce_min_feat_size_func )
            % Synthesizes an apodized grating that radiates desired Gaussian field profile
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
            
            % generate a fiber gaussian mode
            [ obj, field_profile ] = obj.make_gaussian_profile( MFD );
            
            % generate final design, apodized
            if nargin > 3
                obj = obj.generate_final_design_apodized( field_profile, input_wg_type, enforce_min_feat_size_func );
            else
                obj = obj.generate_final_design_apodized( field_profile, input_wg_type );
            end
            
            % save MFD
            obj.synthesized_design.MFD = MFD;
            
        end     % end generate_final_design_apodized_gaussian()
        
        
        function [ obj, chosen_fills, chosen_periods, ...
                    chosen_angles, chosen_scatter_str, ...
                    chosen_ks ] = pick_design_half( obj, input_wg_type )
            % Pick the half of the design space to use depending on whether 
            % input waveguide is full or none
            % 
            % inputs:
            %   input_wg_type
            %       supported options: 'full', 'bottom', 'none' (bottom and
            %       none default to same behavior)
            
            % pick maximum scattering cell and select cells that have larger fill
            % remember fill is sorted in descending order
            [ ~, indx_max_scatter ] = max( obj.sweep_variables.scatter_str_vs_fill );
            
            % choose which half of design space to use, depending on input wg
            switch input_wg_type
                case 'full'
                    % waveguide is full and grating will be etched, so start with smallest gaps (largest fills)
                    chosen_fills        = obj.sweep_variables.fill_ratios_to_sweep( 1:indx_max_scatter );
                    chosen_periods      = obj.sweep_variables.periods_vs_fill( 1:indx_max_scatter );
                    chosen_angles       = obj.sweep_variables.angles_vs_fill( 1:indx_max_scatter );
                    chosen_scatter_str  = obj.sweep_variables.scatter_str_vs_fill( 1:indx_max_scatter );
                    chosen_ks           = obj.sweep_variables.k_vs_fill( 1:indx_max_scatter );
                case {'bottom', 'none'}
                    % grating will be deposited on top of grating, so start with smallest teeth sizes (smallest fills)
                    chosen_fills        = obj.sweep_variables.fill_ratios_to_sweep( indx_max_scatter:end );
                    chosen_periods      = obj.sweep_variables.periods_vs_fill( indx_max_scatter:end );
                    chosen_angles       = obj.sweep_variables.angles_vs_fill( indx_max_scatter:end );
                    chosen_scatter_str  = obj.sweep_variables.scatter_str_vs_fill( indx_max_scatter:end );
                    chosen_ks           = obj.sweep_variables.k_vs_fill( indx_max_scatter:end );
            end
            
        end     % end pick_design_half()
        
        
        function [ obj, chosen_fills, chosen_periods, ... 
                chosen_angles, chosen_scatter_str, chosen_ks ]...
                = enforce_min_feat_size( obj, enforce_min_feat_size_func, ... 
                                            chosen_fills, ... 
                                            chosen_periods, ... 
                                            chosen_angles, ... 
                                            chosen_scatter_str, ... 
                                            chosen_ks )
            % Enforces minimum feature size
            %
            % Inputs:
            %   enforce_min_feat_size_func
            %       function handle
            
            indices_to_keep = [];
            for ii = 1:length( chosen_periods )
                if enforce_min_feat_size_func( chosen_periods(ii), chosen_fills(ii) ) == true
                    indices_to_keep(end+1) = ii;
                end
            end

            chosen_fills        = chosen_fills(indices_to_keep);
            chosen_periods      = chosen_periods(indices_to_keep);
            chosen_angles       = chosen_angles(indices_to_keep);
            chosen_scatter_str  = chosen_scatter_str(indices_to_keep);
            chosen_ks           = chosen_ks(indices_to_keep);                     
            
        end     % end enforce_min_feat_size()
      
        
        function [ obj, field_profile ] = make_gaussian_profile( obj, MFD )
            % Makes gaussian field profile for generate final designs
            % generate a fiber gaussian mode
            w0          = MFD/2;                                                       
            zvec        = 0;                                                            % this is unused
            d0          = 0;                                                            % take slice at waist
            % generate x coordinates for the gaussian mode
            % must be large enough to fit mode
%             xvec_fib        = 0 : obj.discretization : MFD*4.5 - obj.discretization;
            xvec_fib        = 0 : obj.discretization : MFD*6 - obj.discretization;
            xvec_fib        = xvec_fib - xvec_fib(round(end/2));                                % shift origin over to middle
            [obj, field_profile]  = obj.fiber_mode_gaussian(  w0, zvec, xvec_fib, obj.optimal_angle, d0, obj.background_index );
            field_profile         = field_profile./max( abs(field_profile) );
            
        end     % end make_gaussian_profile()
        
        
        function [ obj, overlap_12 ] = xcorr_normalized( obj, field1, field2 )
            % very simple function that just computes normalized cross correlation of two vectors, 
            % field1 and field2
            
            % transpose fields so they are nx1
            field1 = reshape( field1, [ length(field1), 1 ] );
            field2 = reshape( field2, [ length(field2), 1 ] );
            
            % first normalize both fields
            norm_factor_1   = field1' * field1;
            norm_factor_2   = field2' * field2;
            
            % fft may require padding, get length of the longer of the 2 arrays
            array_len = max( [ length(field1), length(field2) ] );
            
            % calculate field1 and field2 xcorr
            overlap_12 = ifftshift( ifft( fft( fftshift( field1 ), array_len ) .* conj( fft( fftshift( field2 ), array_len ) ) ) );
%             overlap_12 = overlap_12./sqrt( ( norm_factor_1 .* norm_factor_2 ) );
            overlap_12 = ( abs(overlap_12).^2 )./( norm_factor_1 .* norm_factor_2 );
            
        end     % end function xcorr_normalized()
        
        
        function [obj, field_shape_prediction] = predict_field_shape( obj, xvec, scatter_str_vs_x )
            % Calculates predicted field amplitude profile by integrating the scattering strength
            %
            % Inputs:
            %   xvec
            %       type: double, array
            %       desc: x vector coordinates 
            %   scatter_str_vs_x
            %       type: double, array
            %       desc: alpha at each x vector coordinate
            %
            % Outputs:
            %   field_shape_prediction
            %       type: double, array
            %       desc: predicted field amplitude profile vs. xvec
            
            % integrate the picked scatter strengths to get the predicted
            % field shape
            field_shape_prediction = zeros( size(xvec) );
            A0 = 1;
            for ii = 1:length(xvec)

                % update radiation
                field_shape_prediction(ii) = A0.*(1 - exp( -scatter_str_vs_x(ii) * obj.discretization ) );

                % update guided power
                A0 = A0.*exp( -scatter_str_vs_x(ii) * obj.discretization );

            end
                
        end     % end predict_field_shape()
        
        function [obj, field_shape_prediction] = predict_field_shape_v2( obj, xvec, scatter_str_vs_x, imag_k_vs_x )
            % Calculates predicted field amplitude profile by integrating the scattering strength
            % updated to use imag(k) to update loss
            %
            % Inputs:
            %   xvec
            %       type: double, array
            %       desc: x vector coordinates 
            %   scatter_str_vs_x
            %       type: double, array
            %       desc: alpha at each x vector coordinate
            %
            % Outputs:
            %   field_shape_prediction
            %       type: double, array
            %       desc: predicted field amplitude profile vs. xvec
            
            % integrate the picked scatter strengths to get the predicted
            % field shape
            field_shape_prediction = zeros( size(xvec) );
            A0 = 1;
            for ii = 1:length(xvec)

                % update radiation
                field_shape_prediction(ii) = A0.*(1 - exp( -scatter_str_vs_x(ii) * obj.discretization ) );

                % update guided power
                A0 = A0.*exp( -imag_k_vs_x(ii) * obj.discretization );

            end
                
        end     % end predict_field_shape()
        
        
        function [ obj, xvec, alpha_des, field_profile ] = calculate_desired_scattering( obj, field_profile )
            % Calculates desired scattering profile for an arbitrary field profile
            %
            % Right now, this function also makes the xvec for you
            %
            % Inputs:
            %   field_profile
            %       type:
            %       desc:
            %
            % Outputs:
            %   xvec
            %   alpha_des
           
            % normalize the field to intensity
            field_profile = field_profile/sqrt( obj.discretization * sum( abs( field_profile ).^2 ) );
            
            % generate x coordinates for the field profile, is this the right place to do this?
            xvec = obj.discretization .* ( 0:length(field_profile)-1 );
            xvec = xvec - xvec(round(end/2));                                % shift origin over to middle
                                              
            % calculate desired scattering strength vs. x
            integral_f      = cumsum( abs(field_profile).^2 ) * obj.discretization;
%             alpha_des       = (1/2)*( abs(u).^2 ) ./ ( 1 + 1e-9 - integral_f );             % in units 1/m
            alpha_des       = (1/2)*( abs(field_profile).^2 ) ./ ( 1 - integral_f );             % in units 1/units
%             alpha_des       = alpha_des * obj.units.scale;                                  % in units 1/units
% 
            % limit the xvec and alpha to certain radiated output threshold
            threshold   = 1 - 1e-30;
            xvec        = xvec( integral_f <= threshold );
            alpha_des   = alpha_des( integral_f <= threshold );
            field_profile = field_profile( integral_f <= threshold );
            
            % save the xvec, alpha, and field as variables of final
            % syntheiszed design
            obj.synthesized_design.xvec         = xvec;
            obj.synthesized_design.alpha_des    = alpha_des;
            obj.synthesized_design.field_profile = field_profile;
            
        end     % end function calculate_desired_scattering()
        
        function [ obj, synthesized_design ] = pick_final_datapoints( obj, ...
                                      xvec, alpha_des, ...
                                      chosen_fills, ...
                                      chosen_periods, ...
                                      chosen_angles, ...
                                      chosen_scatter_strs, ...
                                      chosen_ks, ...
                                      start_alpha_des )
            % Picks the final datapoints (cells) that make up the grating
            %
            % Inputs:
            %     xvec
            %         type: double, array
            %         desc: coordinates of desired field profile
            %     alpha_des
            %         type: double, array
            %         desc: desired scattering strength vs. x
            %     chosen_fills
            %         type: double, array
            %         desc: swept fills
            %     chosen_periods
            %         type: double, array
            %         desc: period vs. fill from generated design space
            %     chosen_angles
            %         type: double, array
            %         desc: angle vs fill from generated design space
            %     chosen_scatter_strs
            %         type: double, array
            %         desc: scattering strength (alpha) from generated design space
            %     chosen_ks
            %         type: double, array
            %         desc: k from generated design space
            %     start_alpha_des
            %         type: double, scalar
            %         desc: OPTIONAL starting alpha
            %
            % Outputs:
            %   synthesized_design:
            %       type: struct
            %       desc: a struct that holds the synthesized cells
            
            % default to start from weakest cell
            if ~exist( 'start_alpha_des', 'var' )
                start_alpha_des = min(chosen_scatter_strs);
            end

            % now match these data points to the desired alpha
            % starting point
            [~, indx_max_alpha] = max( alpha_des );
            [~, indx_x]         = min(abs( alpha_des(1:indx_max_alpha) - start_alpha_des ) );
            cur_x               = xvec(indx_x);
            
            % final synthesized variables
            synthesized_design.fill         = [];
            synthesized_design.period       = [];
            synthesized_design.angles       = [];
            synthesized_design.scatter_str  = [];
            synthesized_design.k            = [];
            synthesized_design.GC           = {};
            synthesized_design.des_scatter  = [];
            
            % flag for switching to using max scattering strength
            saturate_scatter_str_to_max = false;
 
            ii = 1;
            while cur_x < xvec(end)
                % build grating one cell at a time
                
                % pick design with scattering strength closest to desired
                % alpha
                des_scatter = alpha_des(indx_x);                            % desired alpha
                if des_scatter  > max( chosen_scatter_strs )
                    % desired scattering strength too high, gotta saturate
                    saturate_scatter_str_to_max = true;
                end
                if ~saturate_scatter_str_to_max
                    [~, indx_closest_scatter]   = min( abs(chosen_scatter_strs - des_scatter) );          % index of closest scatter design 
                else
                    [~, indx_closest_scatter]   = max( chosen_scatter_strs );                             % saturate to max
                end
                
                % save parameters
                synthesized_design.fill(ii)         = chosen_fills( indx_closest_scatter );
                synthesized_design.period(ii)       = chosen_periods( indx_closest_scatter );
                synthesized_design.angles(ii)       = chosen_angles( indx_closest_scatter );
                synthesized_design.scatter_str(ii)  = chosen_scatter_strs( indx_closest_scatter );
                synthesized_design.k(ii)            = chosen_ks( indx_closest_scatter );
                synthesized_design.des_scatter(ii)  = des_scatter;
                
                % save a grating cell object for each cell
                synthesized_design.GC{ii} = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               synthesized_design.period(ii), ...
                                               synthesized_design.fill(ii) );
                
                % move onto next
                cur_x       = cur_x + synthesized_design.period(ii);
                [~, indx_x] = min( abs(xvec - cur_x) );
                cur_x       = xvec( indx_x );
                ii          = ii + 1;
                
            end     % end for ii = 1:ncells
            
        end     % end pick_final_datapoints()
        
        
        function obj = build_final_index( obj )
            % build final index distribution
            % to be run after a design has been synthesized
            %
            % currently directly modifies obj.synthesized_design
        
            obj.synthesized_design.N = [];
            for ii = 1:length(obj.synthesized_design.GC)

                GC                          = obj.synthesized_design.GC{ii};
                obj.synthesized_design.N    = [ obj.synthesized_design.N, GC.N ];

            end
            
            % coordinates of index distribution
            obj.synthesized_design.x_coords = obj.discretization*( 0:1:( size(obj.synthesized_design.N,2)-1 ) );
            obj.synthesized_design.y_coords = obj.discretization*( 0:1:( size(obj.synthesized_design.N,1)-1 ) );
            
        end     % end build_final_index()
        
        
        function [ obj, best_alpha_power ] = optimize_start_alpha( obj, xvec, alpha_des,  ...
                                                chosen_fills, ...
                                                chosen_periods, ...
                                                chosen_angles, ...
                                                chosen_scatter_strs, ...
                                                chosen_ks, ...
                                                desired_field )
            % Optimizes the starting alpha, based on predicted field overlap
            %
            % Inputs:
            %     desired_field
            %         type: double, array
            %         desc: desired field vs. xvec
           
            % using matlab's fminsearch
            f = @(alpha_power)( 1 - obj.predict_overlap_for_optimization( ...
                                              xvec, alpha_des,  ...
                                                chosen_fills, ...
                                                chosen_periods, ...
                                                chosen_angles, ...
                                                chosen_scatter_strs, ...
                                                chosen_ks, ...
                                                desired_field, alpha_power ) );
                  
            best_alpha_power = fminsearch( f, log10(min( chosen_scatter_strs )) );
                                                   
        end     % end function optimize_start_alpha()
        
        
        function [ max_overlap, field_shape_prediction ] = predict_overlap_for_optimization( ...
                                              obj, xvec, alpha_des, ...
                                                chosen_fills, ...
                                                chosen_periods, ...
                                                chosen_angles, ...
                                                chosen_scatter_strs, ...
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
            [ obj, synthesized_design ] = obj.pick_final_datapoints( xvec, alpha_des, ...
                                              chosen_fills, ...
                                              chosen_periods, ...
                                              chosen_angles, ...
                                              chosen_scatter_strs, ...
                                              chosen_ks, ...
                                              10.^(start_alpha_power));

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
            [obj, field_shape_prediction]   = obj.predict_field_shape( xvec, scatter_str_vs_x );
                
            % calculate overlap with desired field (simple xcorr)
            [ ~, overlap ]  = obj.xcorr_normalized( desired_field, field_shape_prediction );
            max_overlap     = max( abs(overlap) );
            
        end     % end function predict_overlap_for_optimization()
        
    end     % End methods section
    
end     % end class definition



% -------------------------------------------------------------------------
% Begin auxiliary non-class methods
% -------------------------------------------------------------------------

function GC = makeGratingCell( synth_obj, period, fill_top, fill_bot, offset_ratio )
% currently deprecated
% makes and returns a c_twoLevelGratingCell object
% 
% inputs:
%   synth_obj
%       type: c_synthGrating object AS STRUCT
%       desc: c_synthGrating object AS STRUCT
%   period
%       type: double, scalar
%       desc: period of the grating cell
%   fill_top
%       type: double, scalar
%       desc: ratio of top layer to period
%   fill_bot
%       type: double, scalar
%       desc: ratio of bottom layer to bottom layer
%   offset_ratio
%       type: double, scalar
%       desc: ratio of bottom layer offset to period
%
% outputs:
%   GC
%       type: c_twoLevelGratingCell object
%       desc: two level grating cell object


% set domain 
domain_size     = synth_obj.domain_size;
domain_size(2)  = period;

% make grating cell
GC = c_twoLevelGratingCell( 'discretization', synth_obj.discretization, ...
                            'units', synth_obj.units.name, ...
                            'lambda', synth_obj.lambda, ...
                            'domain_size', domain_size, ...
                            'background_index', synth_obj.background_index );

                        
% wrap offsets to range 0 to 1
offset_ratio = mod( offset_ratio, 1 );
                        
% draw cell
% draw two levels using two level builder function
% the inputs are organized [ top level, bottom level ]
wg_thick        = 100.0 / ( 1e9 * synth_obj.units.scale ) * ones(1,2);      % hardcoded to 100nm
wg_index        = [3.45, 3.45];                                             % hardcoded to 3.45
wg_min_y        = [ domain_size(1)/2, domain_size(1)/2-wg_thick(1) ];
wgs_duty_cycles = [ fill_top, fill_bot ];
wgs_offsets     = [ 0, offset_ratio*period ];
GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, wg_index, ...
                                        wgs_duty_cycles, wgs_offsets );
            
end




















































