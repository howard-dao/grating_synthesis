classdef c_gratingCell < c_bloch_cell
% Class for simulating a grating cell, using FDFD complex-k solver
%
% Authors: bohan zhang
%
% Description:
%
% Notes:
%   - propagation direction coordinate is 'x'. 
%   - transverse (in plane) direction coordinate is 'y'.
%   - transverse (out of plane) direction is 'z'.
%
% Inputs to constructor:
%   Inputs are name-value pairs:
%   'discretization'
%       type: double, scalar or 1x2 vector
%       desc: discretization along x and y, in units of 'units'
%             if scalar, then dx = dy = discretization
%             if vector, then discretization = [ dy dx ]
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
%   'num_cells'
%       type: integer, scalar
%       desc: OPTIONAL: pick the number of grating cells to repeat when
%             drawing/calculating Ez
%
% Example usage:
%
%   % make grating cell
%   GC = c_gratingCell( 'discretization', dxy, ...
%                       'units', units, ...
%                       'lambda', lambda, ...
%                       'domain_size', [ y_domain_size, period ], ...
%                       'background_index', background_index, ...
%                       'numcells', 10 );
    
    properties
        
        % geometry properties
        wg_min_y;       % bottom position of wg
        wg_max_y;       % top position of wg
        
        % radiation parameters
        directivity;        % ratio of power radiated up/power radiated down
        max_angle_up;       % Angle of maximum upwards output radiation
        max_angle_down;     % Angle of maximum downwards output radiation
        P_rad_down;         % saves power radiated down
        P_rad_up;           % saves power radiated up
        P_in;               % saves power entering the wg
        alpha_up;           % fractional radiative loss/unit length upwards
        alpha_down;         % fractional radiative loss/unit length downwards
        Sx;                 % X component of poynting vector (dims y vs. x(2:end-1))
        P_per_x_slice;      % x traveling power per x slice vs. x(2:end-1)
        Sy;                 % Y component of poynting vector (dims y(2:end-1) vs x)
        P_per_y_slice;      % Y traveling (upwards) power per y slice vs. y(2:end-1)
        % some new additions, 8/18/19
        Sy_up           
        Sy_down      
        P_thru      
        alpha_up_from_srad
        % added 12/5/20
        n_top;
        n_bot;
        
    end     % end properties
    
    
    methods
        
        % -----------------------------------------------------------------
        % Constructor
        
        function obj = c_gratingCell( varargin )
            
            % call bloch cell constructor
            obj = obj@c_bloch_cell(varargin{:});
            
        end     % end constructor
        
        % -----------------------------------------------------------------
        % Running simulations
        
        function obj = runSimulation( obj, num_modes, BC, pml_options, k0, guessk, OPTS )
            % Runs new mode solver
            %
            % Description:
            %   Runs complex-k mode solver. Stores the mode with the most
            %   guided power. Calculates up/down power, directivity,
            %   angle of maximum radiation, and scattering strength.
            %
            % Inputs:
            %   num_modes
            %       type: integer
            %       desc: # of modes (max) to simulate
            %   BC
            %       type: integer
            %       desc: 0 for PEC, 1 for PMC
            %   pml_options
            %       type: array, double
            %       desc: 1x4 Array with the following elements:
            %               PML_options(1): PML in y direction (yes=1 or no=0)
            %               PML_options(2): length of PML layer in nm
            %               PML_options(3): strength of PML in the complex plane
            %               PML_options(4): PML polynomial order (1, 2, 3...)
            %   k0
            %       type: scalar, double
            %       desc: free-space wavenumber
            %   guessk
            %       type: scalar, double (can be complex)
            %       desc: guess k value. Works best when closest to desired
            %             mode. In units rad/'units'
            %   OPTS
            %       type: struct
            %       desc: optional options with the following fields
            %           'mode_to_overlap'
            %               type: EITHER a matrix, double, or a scalar, int
            %               desc: mode to overlap (if matrix), or # of mode
            %                       to pick if scalar
            %
            % Sets these properties:
            %   obj.k
            %       units rad/'units'
            %   obj.Phi
            %       dimensions y vs. x (x is dir. of propagation)
            %   obj.E_z
            %       field repeated, i can't remember why or if i use this
            %   obj.directivity
            %       up/down power ratio
            %   calls obj.calc_output_angle()
            %       this is a function, which probably figures out the
            %       output angle
            %   calls obj.calc_scattering_strength()
            
            % default OPTS
            if nargin < 7
                OPTS = struct();
            end
            
            % first run the bloch cell simulation function
            obj = runSimulation@c_bloch_cell( obj, num_modes, BC, pml_options, k0, guessk, OPTS );
            
            % pick which mode to keep
            if isfield( OPTS, 'mode_to_overlap' )
                % pick mode with best overlap
                obj = obj.choose_mode( OPTS.mode_to_overlap );
            else
                % pick mode guided mode
                obj = obj.choose_mode(); 
            end
            
            % stitch together full e field, with the request number of
            % periods
            % i don't even remember what this is for
            % looks like it's used in calc_scattering_strength
            [obj, obj.E_z]  = obj.stitch_E_field( obj.Phi, obj.k, obj.numcells );
            
            % for mode overlapping, stitch together single unit cell of E
            % field, using only real(k)
            [obj, E_z_for_overlap]  = obj.stitch_E_field( obj.Phi, real(obj.k), 1 );
            obj.E_z_for_overlap     = E_z_for_overlap;
            
            % pick slices of field to compute directivity, angle, etc.
            h_pml_d = round( pml_options(2)/obj.dy );                       % size of pml in discretizations
            y_up    = size( obj.E_z, 1 ) - h_pml_d - 1;
            y_down  = h_pml_d+2;
            
            % calculate output angle
            obj = obj.calc_output_angle(k0);
            
            % calculate up/down directivity
            obj             = obj.calc_radiated_power(k0);
            obj.directivity = obj.P_rad_up/obj.P_rad_down;
                    
            % calculate power scattering strength
            % alpha is imaginary k
            obj.alpha_up    = imag(obj.k);
            obj.alpha_down  = imag(obj.k);
            
        end     % end runSimulation
        
        
        function obj = calc_output_angle(obj, k0)
            % Calculates output angle of radiation
            % To be run only AFTER mode solver has been run

            pml_size        = round(obj.sim_opts.pml_options(2)/obj.dy);
            
            % calculate angle from phase matching
            kx      = real(obj.k) - 2*pi/obj.domain_size(2);
            n_top   = obj.N( end-pml_size-1, 1 );
            n_bot   = obj.N( pml_size+2, 1 );
            k0_cl_up   = n_top * k0;
            k0_cl_down = n_bot * k0;
            obj.max_angle_up    = asin( kx./k0_cl_up ) * 180/pi;
            obj.max_angle_down  = asin( kx./k0_cl_down ) * 180/pi;
            
        end     % end function calc_output_angle()
        
        
        function obj = calc_radiated_power(obj, k0)
            % Calculates the radiated power in both the upwards and
            % downwards directions
            % Only run AFTER mode solver has been run, since this function
            % depends on E_z and k
            %
            % Inputs:
            %   y_up    - location of top slice of E field to take
            %   y_down  - location of bot slice of E field to take
            %       y_up and y_down are deprecated.
            %
            % Uses these object properties:
            %   E_z
            %   k
            %
            % Sets these object properties:
            %   P_rad_down      - saves power radiated down
            %   P_rad_up        - saves power radiated up
                       
            % calculate power distributions for the chosen mode
            [obj, Sx, Sy, P_per_x_slice, P_per_y_slice] = obj.calc_power_distribution_ncells(1, k0);
            
            % take poynting vector at vertical boundaries
            pml_size = round(obj.sim_opts.pml_options(2)/obj.dy);
            Sy_up    = Sy( end-pml_size, : ); % dimensions Sy vs. x
            Sy_down  = Sy( pml_size+1, : );
            
            % save to obj
            obj.Sx              = Sx;
            obj.P_per_x_slice 	= P_per_x_slice;
            obj.Sy              = Sy;
            obj.Sy_up           = Sy_up;
            obj.Sy_down         = Sy_down;
            obj.P_per_y_slice   = P_per_y_slice;
            
            % radiated power is integral of poynting vector
            obj.P_rad_up    = sum( abs( Sy_up(:) ) ) * obj.dx ;
            obj.P_rad_down  = sum( abs( Sy_down(:) ) ) * obj.dx ;
            obj.P_in        = sum( abs( Sx( :, 1 ) ) ) * obj.dy ;
            obj.P_thru      = sum( abs( Sx( :, end ) ) ) * obj.dy ;
     
        end     % end function calc_radiated_power()
        
        
        function obj = calc_scattering_strength(obj)
            % deprecated
            % Calculates scattering strength in terms of decay rate alpha
            % To be run only AFTER mode solver AND calc_radiated_power have
            % been run
            
            % define constants
            mu0     = 4*pi*1e-7;                % units of H/m
            c       = 3e8;                      % units of m/s
            omega0 	= 2*pi*c/obj.lambda;        % units of rad*(m/s)/units
            
            % grab field
           
            % how about H y in, on FIRST step, using periodicity
            
            % stitch together e field, including the phase
            % on the -1, 0, and 1 steps
            E_z_neg1        = obj.Phi( :, end ) .* exp( 1i * obj.k * -obj.dx );
            E_z_0           = obj.Phi( :, 1 );
            E_z_plus1       = obj.Phi( :, 2 ) .* exp( 1i * obj.k * obj.dx );  
            
            H_y_in = (1i/(omega0*mu0)) * ( E_z_plus1 - E_z_neg1 )/(2*obj.dx);   % dy, term staggered 1
            
            % power in (at left edge)
            Sx_in   = real( -1 * conj(H_y_in(:)) .* E_z_0 );              % Sx = real( -Ez Hy* )
            P_in    = sum( Sx_in(:) )*obj.dy;
            
            
            % DEBUG save Sx_in
            obj.debug.Sx_in = Sx_in;

            % save input power
            obj.P_in = P_in;
            
            
            
            % calculate alpha efficiency (in 1/units)
            period          = obj.domain_size(2);
            
            % an old version
%             obj.alpha_up    = ( -1/( 2*period ) ) * log( 1 - obj.P_rad_up/P_in );
%             obj.alpha_down  = ( -1/( 2*period ) ) * log( 1 - obj.P_rad_down/P_in );

            % another old version
            % this one calcualtes it using decay length of intensity of radiated beam, but
            % only across a single period.
%             obj.alpha_up    = ( -1/( 2*period ) ) * log( obj.Sy_up(end)/obj.Sy_up(1) );
%             obj.alpha_down  = ( -1/( 2*period ) ) * log( obj.Sy_down(end)/obj.Sy_down(1) );
%             % this one calculates it using ratio of radiated power to input
%             % power minus the other losses
%             Prad_in_up = obj.P_in - obj.P_rad_down - obj.P_thru;
%             alpha_up_from_Prad = ( -1/( 2*period ) ) * log( 1 - obj.P_rad_up/Prad_in_up );

            % newer ver, 6/7/2020 that first predicts how many periods to
            % look at
            
%             % first calculate decay length, so power decays to e^(-4) =
%             % 1.8%
%             decaylen = 2/imag(obj.k);
%             
%             % now pick num of periods
%             n_periods = round( decaylen/period );
%             if isinf(n_periods)
%                 n_periods = 1;
%             end
%             [~, Sx, Sy, P_per_x_slice, P_per_y_slice] = obj.calc_power_distribution_ncells(n_periods);
%             pml_size        = round(obj.sim_opts.pml_options(2)/obj.dy);
%             obj.alpha_up    = ( -1/( 2*n_periods*period ) ) * log( Sy(end - pml_size, end)/Sy(end - pml_size, 1) );
%             obj.alpha_down  = ( -1/( 2*n_periods*period ) ) * log( Sy(pml_size+1, end)/Sy(pml_size+1, 1) );
                     
            
%             % new ver, 6/24/2020 that uses total radiated power to get
%             % alpha
%             % first calculate decay length, so power decays to e^(-4) =
%             % 1.8%
%             decaylen = 2/imag(obj.k);
%             
%             % now pick num of periods
%             n_periods = ceil( decaylen/period );
%             if isinf(n_periods)
%                 n_periods = 1;
%             end
%             
%             % if too many periods to fit into memory, slim it down
%             max_size = 1e5;
%             if n_periods*period/obj.dx > max_size
%                 n_periods = floor( max_size*obj.dx/period );
%             end
%             
%             [~, Sx, Sy, P_per_x_slice, P_per_y_slice] = obj.calc_power_distribution_ncells(n_periods);
%             pml_size        = round(obj.sim_opts.pml_options(2)/obj.dy);
%             
%             P_rad_up    = abs(sum( Sy(end - pml_size,:) ) * obj.dx);
%             P_rad_down  = abs(sum( Sy(pml_size+1,:) ) * obj.dx); 
          
%             obj.alpha_up    = ( -1/( 2*period*n_periods ) ) * log( 1 - (P_rad_up+P_rad_down)/P_in );
%             obj.alpha_down  = obj.alpha_up;
            

            
%             % DEBUG let's calculate, numerically, the loss coefficient from
%             % input to output and compare with imag(k)
%             % H x out, at grid point end-1
%             H_x_out1 = (-1i/(omega0*mu0)) * ( E_z( 3:2:end, end-1 ) - E_z( 1:2:end-2, end-1 ) )/(2*obj.dy);   % dy, term staggered 1
%             H_x_out2 = (-1i/(omega0*mu0)) * ( E_z( 4:2:end, end-1 ) - E_z( 2:2:end-2, end-1 ) )/(2*obj.dy);   % dy, term staggered 2
%             % interleave two arrays
%             H_x_out                                 = zeros( size( [H_x_out1, H_x_out2] ) );
%             H_x_out( ( 1:length(H_x_out1) )*2 - 1)  = H_x_out1;
%             H_x_out( ( 1:length(H_x_out2) )*2 )     = H_x_out2;
%             
%             % H y out (on the same grid as Hx out)
%             H_y_out = (1i/(omega0*mu0)) * ( E_z( 2:end-1, end ) - E_z( 2:end-1, end-2 ) )/(2*obj.dx);   % dy, term staggered 1
%             
%             % power in (at left edge)
%             Sy_out = real( conj(H_x_out(:)) .* E_z( 2:end-1, end-1 ) );         % Sy = real( Ez Hx* )
%             Sx_out = real( -1 * conj(H_y_out(:)) .* E_z( 2:end-1, end-1 ) );    % Sx = real( -Ez Hy* )
%             P_out  = sum( sqrt( Sy_out.^2 + Sx_out.^2 ) )*obj.dy;                 % using both Sx and Sy compmonents
%             
%             % numerical calculation of loss coeff
%             alpha_tot = ( -1/( 2*obj.numcells*period ) ) * log( P_out/P_in );
            
        end     % end function calc_scattering_strength()
        
        function [obj, Sx, Sy, P_per_x_slice, P_per_y_slice] = ...
                        calc_power_distribution_ncells(obj, n_periods, k0, mode_num)
            % Calculates Sx and Sy, power distribution for arbitrary number of cells
            %
            % Pre-requisite for running calc_radiated_power
            %
            % currently also saves H fields...
            %
            % inputs:
            %   n_periods
            %       type: scalar, int
            %       desc: number of cells/periods to calculate Sx and Sy for
            %   mode_num
            %       type: scalar, int
            %       desc: OPTIONAL, mode number to calc Sx and Sy for
            %
            % outputs:
            %   Sx
            %       type: matrix, double
            %       desc: x component (in dir. of waveguide prop) of poynting vector
            %             dimensions y vs. x(2:end-1)
            %   Sy
            %       type: matrix, double
            %       desc: y component (transverse dir, pos = up) of poynting vector
            %             dimensions y(2:end-1) vs. x
            %   P_per_x_slice
            %       type: vector, double
            %       desc: integral of Sx aka power propagating through each
            %             x slice of the domain
            %             dimensions power vs. x(2:end-1)
            %   P_per_y_slice
            %       type: vector, double
            %       desc: integral of Sy aka power propagating through each
            %             y slice of the domain
            %             dimensions power vs. y(2:end-1)
            
            if nargin < 4
                choose_mode = false;
            else
                choose_mode = true;
            end
            
            % coordinates
            xcoords = 0 : obj.dx : round((n_periods*obj.domain_size(2) - obj.dx));
            
            % Choose field to use
            if ~choose_mode
                % use the pre-chosen mode

                % Stich field and phase together
                phase_manycells     = repmat( exp( 1i * obj.k * xcoords ), size( obj.Phi, 1 ), 1 );
                Phi                 = repmat( obj.Phi, 1, n_periods );
                Ez_manycells        = Phi .* phase_manycells;
                
            else
                % choose which mode to use
                
                k_chosen    = obj.k_vs_mode( mode_num );
                phi_chosen  = obj.Phi_vs_mode( mode_num );
                % Stich field and phase together
                phase_manycells   = repmat( exp( 1i * k_chosen * xcoords ), size( phi_chosen, 1 ), 1 );
                Phi               = repmat( phi_chosen, 1, n_periods );
                Ez_manycells      = Phi .* phase_manycells;
                
            end
       
            % Calc power vs. x
            
            % add to Ez_onecell Ez( -dx ) and Ez( end + dx ) to calculate H
            E_z_neg1        = Phi( :, end ) .* exp( 1i * obj.k * (xcoords(1)-obj.dx) );
            E_z_plus        = Phi( :, 1 ) .* exp( 1i * obj.k * (xcoords(end)+obj.dx) );
            Ez_wbounds      = [ E_z_neg1, Ez_manycells, E_z_plus ];
            
            % H y in , on second to end-1 steps, using entire E_z
            % dimensions are y (transverse) vs. x
            H_y = (1i/(k0 * obj.c * obj.mu0)) * ( Ez_wbounds( :, 3:end ) - Ez_wbounds( :, 1:end-2 ) )/(2*obj.dx);   % dx, on same grid as Ez_onecell
            
            % poynting vector in, dimensions are y vs x(2:end-1)
            Sx = real( -1 * conj( H_y ) .* Ez_manycells ); % Sx = real( -Ez Hy* )

            % power per x slice
            P_per_x_slice = sum( Sx, 1 )*obj.dy;
            
            % Calculate power radiated vs. y
            
            % calculate H_x
            % dimensions H_x vs. y vs x
            H_x  = 1/(1i * k0 * obj.c * obj.mu0) .* ( Ez_manycells( 3:end,:) - Ez_manycells( 1:end-2,:) )/(2*obj.dy);

            % power flowing across y
            Sy                  = real( conj( H_x ) .* Ez_manycells( 2:end-1,: ) );       % Sy = real( Ez Hx* )
            P_per_y_slice       = sum( Sy, 2 )*obj.dx;                                  % using Sy compmonent
                      
        end     % end function calc_power_distribution_onecell()
        
        
        function obj = choose_mode( obj, mode )
            % Function that chooses which mode becomes the accepted mode
            %
            % Inputs:
            %   mode
            %       type: matrix, double OR scalar integer
            %       desc: OPTIONAL. If passed in as argument, this function
            %             will choose the mode with the closest overlap.
            %             Otherwise, the function will choose the mode with
            %             the mode guided power
            %             if a scalar integer, selects that mode as the
            %             primary mode
            %
            % Sets these properties:
            %   obj.k
            %       prop. k of chosen mode
            %   obj.Phi
            %       field envelope of chosen mode
            %   obj.chosen_mode_num
            %       index of chosen mode (chosen mode k = obj.k_vs_mode(
            %       obj.chosen_mode_num))
            
            if nargin < 2
                % sort on guided power
                
                guided_power  = zeros( size(obj.k_vs_mode) );
                total_power   = zeros( size(obj.k_vs_mode) );

                % check to see if waveguide boundaries have been set yet
                if isempty(obj.wg_min_y) || isempty(obj.wg_max_y)
                    error(['Waveguide boundaries have not been set yet. You must set the waveguide boundaries by either calling' ...
                            ' "twoLevelBuilder()" or by setting the "wg_min_y" AND "wg_max_y" object properties yourself.']);
                end

                y_bot   = obj.wg_min_y;
                y_top   = obj.wg_max_y;
                y       = obj.y_coords;

                % Calculated guided power and fix backwards modes if option is
                % selected
                for ii = 1:length(obj.k_vs_mode)
                    % for each mode

                    % grab guided portion of field
                    phi_guided  = obj.Phi_vs_mode( y >= y_bot & y <= y_top, :, ii );
                    cur_phi     = obj.Phi_vs_mode( :, :, ii );

                    % sum area of field
                    guided_power(ii)    = sum( abs( phi_guided(:) ).^2 );
                    total_power(ii)     = sum( abs( cur_phi(:).^2 ) );
                    
                    % if the imag(k) is negative then discard it from the
                    % options by setting guided power to -1
                    if imag(obj.k_vs_mode(ii)) < 0
                        guided_power(ii) = -1;
                    end

                end     % end for

                % DEBUG storing the guided power
                obj.debug.guided_power          = guided_power;
                obj.debug.guided_power_ratio    = guided_power./total_power;

                % keep mode with LEAST unguided power OR MOST guided power
                [~, indx_k]         = max( guided_power./total_power );        % most guided
                obj.k               = obj.k_vs_mode(indx_k);
                obj.Phi             = obj.Phi_vs_mode(:,:,indx_k);
                obj.chosen_mode_num = indx_k;
                
            else
                if isscalar( mode )
                    % pick which mode to use as the primary mode
                    obj.k               = obj.k_vs_mode(mode);
                    obj.Phi             = obj.Phi_vs_mode(:,:,mode);
                    obj.chosen_mode_num = mode;
                    
                else
                    % sort on mode overlap
                    obj = choose_mode@c_bloch_cell( obj, mode );
                end
                
            end
            
        end     % end function choose_mode()
        
        
        function obj = calc_directivity_single_angle( obj )
            % Calculate directivity based on a single angle
            % Requires that simulation has already been run
            %
            % this function is very much a WIP, not used for synthesis
            %
            % inputs:
            
            if nargin < 2
                choose_mode = false;
            else
                choose_mode = true;
            end
            
            % define constants
            mu0     = 4*pi*1e-7;                % units of H/m
            c       = 3e8;                      % units of m/s
            omega0 	= 2*pi*c/obj.lambda;        % units of rad m s^-1 units^-1
       
            % Stich field and phase together
            Ez_onecell      = obj.Phi .* repmat( exp( 1i * obj.k * obj.x_coords ), size( obj.Phi, 1 ), 1 );
            
            % calculate H_x
            % dimensions H_x vs. y vs x
            H_x  = 1/(1i*omega0*mu0) .* ( Ez_onecell( 3:end,:) - Ez_onecell( 1:end-2,:) )/(2*obj.dy);

            % pick monitor planes
            [~, indx_max_up]    = max( obj.P_per_y_slice );
            [~, indx_max_down]  = min( obj.P_per_y_slice );
         
        end     % end function calc_directivity_single_angle()
        
        
    end     % end methods
    
    
end

