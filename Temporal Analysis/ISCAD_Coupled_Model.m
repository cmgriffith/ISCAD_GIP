classdef ISCAD_Coupled_Model < matlab.System
    % ISCAD_Coupled_Model_Fast
    % High speed Coupled Circuit Model using pre calculated LUTs.
    % Automatically aligns grid points to stator/rotor slotting to capture peaks.
    
    % Public Properties
    properties (Nontunable)
        mparams = struct();
        % LUT Resolution Multiplier
        % 1 = Points exactly at slot alignments. 10 = 10 points between alignments.
        Resolution_Factor = 10; 
        f_sw = 20e3; % Switching Frequency [Hz]
        Qs = 60; % Number of Stator Slots
        Qr = 50; % Number of Rotor Slots
        p = 2; % Number of Poles

        rho_s = 0; % Resistivity of Conductor
        h_slot = 0; % Slot Height 
        ag = 0; 
        l_stack = 0; % Stack Length
        w_slot = 0; % Slot Width
        r = 0; % Bore Radius
        h_so = 0; % Slot Opening Height
        w_so = 0; % Slot Opening Width

        Ts = 1e-6;  % Step time
    end

    properties (Access = private)
        % Derived parameters for internal use
        Rs_DC, Rr_DC
        L_ss, L_ls, L_lr
    
        InvL_LUT % 3D Matrix: [110 x 110 x Steps] Stores inv(L)
        G_LUT % 3D Matrix: [110 x 110 x Steps] Stores dL/dTheta
        R_Matrix  % Constant resistance matrix
        
        d_theta   % Step size of the LUT (radians)
        Num_Steps % Number of steps in 360 degrees
    end

    properties(DiscreteState)
        Currents 
    end

    methods(Access = protected)
        
        function setupImpl(obj)
            % 0. Fixed Calculations
            % Constants
            mu_0 = 4*pi*1e-7;
            
            % Calculate Resistances (DC)
            % Stator
            % Note: I think w_slot is [inner, outer]. Using index 1.
            obj.Rs_DC = obj.rho_s * obj.l_stack / (obj.h_slot * obj.w_slot(1));
            
            % Rotor (Assuming similar area scaling as analytical script)
            A_rotorbar = obj.h_slot * obj.w_slot(1) * obj.Qs / obj.Qr;
            obj.Rr_DC = obj.rho_s * obj.l_stack / A_rotorbar;
            D_rotorbar = sqrt(4*A_rotorbar/(pi)); % Equivalent diamater of each rotor bar
            
            % Calculate Inductances
            % Self Inductance (Approximate Peak)
            obj.L_ss = (pi/2) * mu_0 * (obj.r * obj.l_stack / obj.ag);
            
            % Leakage Inductance (Only using fixed switching frequency)
            zeta_Ls = sqrt(pi*mu_0*obj.f_sw/obj.rho_s)*obj.h_slot;
            K_Ls = (3/(2*zeta_Ls))*(sinh(2*zeta_Ls)-sin(2*zeta_Ls))/(cosh(2*zeta_Ls)-cos(2*zeta_Ls));
            zeta_Lr = sqrt(pi*mu_0*obj.f_sw/obj.rho_s)*D_rotorbar;
            K_Lr = (3/(2*zeta_Lr))*(sinh(2*zeta_Lr)-sin(2*zeta_Lr))/(cosh(2*zeta_Lr)-cos(2*zeta_Lr));
            
            lambda_slot = (obj.h_slot / (3*obj.w_slot(1))) + (obj.h_so / obj.w_so);
            obj.L_ls = mu_0 * obj.l_stack * lambda_slot;
            obj.L_lr = obj.L_ls; % Assumption for rotor
            
            % Local vars for LUT generation
            L_ss = obj.L_ss;
            L_ls = obj.L_ls;
            L_lr = obj.L_lr;
            
            % 1. Determine Optimal Discretisation
            % Calculate slot pitches in degrees
            deg_s = 360 / obj.Qs;
            deg_r = 360 / obj.Qr;
            
            % Find the greatest common divisor (GCD) of the pitches
            % This ensures the grid hits exactly where stator and rotor align.
            % (Using a small tolerance for floating point GCD)
            common_pitch = gcd_float(deg_s, deg_r); 
            
            % Define step size
            step_deg = common_pitch / obj.Resolution_Factor;
            obj.Num_Steps = round(360 / step_deg);
            obj.d_theta = deg2rad(step_deg);
            
            % Initialise State
            obj.Currents = zeros(obj.Qs + obj.Qr, 1);
            
            % 3. Pre Calculate LUTs
            fprintf('Pre-calculating %d matrices', obj.Num_Steps);
            
            % Pre-allocate for speed
            N_sys = obj.Qs + obj.Qr;
            obj.InvL_LUT = zeros(N_sys, N_sys, obj.Num_Steps);
            obj.G_LUT    = zeros(N_sys, N_sys, obj.Num_Steps);
            
            % Stator/Rotor angle vectors (0 position)
            as_base = (0:obj.Qs-1)' * (2*pi/obj.Qs);
            ar_base = (0:obj.Qr-1)' * (2*pi/obj.Qr);
            
            for k = 1:obj.Num_Steps
                theta = (k-1) * obj.d_theta;
                
                % Build L(theta)
                % Stator-Stator (Fixed)
                Lss = L_ss * cos(as_base - as_base') + eye(obj.Qs)*L_ls;
                
                % Rotor-Rotor (Fixed in Rotor Frame)
                Lrr = L_ss * cos(ar_base - ar_base') + eye(obj.Qr)*L_lr;
                
                % Mutual (Varies with Theta)
                % Angle diff = StatorAngle - (RotorAngle + Theta)
                ar_now = ar_base + theta;
                angle_sr = as_base - ar_now';
                Msr = L_ss * cos(angle_sr);
                
                % Assemble L
                L_mat = [Lss, Msr; Msr', Lrr];
                
                % Build G(theta) = dL/dTheta
                % Only Msr changes with theta. d/dTheta(cos(a - theta)) = sin(a - theta)
                dMsr = L_ss * sin(angle_sr);
                G_mat = [zeros(obj.Qs), dMsr; dMsr', zeros(obj.Qr)];
                
                % Store in LUT
                % Store inverse L to save computing inv() at runtime
                obj.InvL_LUT(:,:,k) = inv(L_mat);
                obj.G_LUT(:,:,k)    = G_mat;
            end
            fprintf('Done.\n');
        end
        
        function [I_stator, Torque] = stepImpl(obj, V_stator, Theta, Omega, f_elec)
            Qs = obj.Qs;
            Qr = obj.Qr;
            
            % 1. Calculate dynamic resistance using fundimental frequency
            if f_elec < 1 % Avoid division by zero if rotor is at standstill
                K_Rs = 1; % No additional AC resistance at zero speed
                K_Rr = 1;
            else
                % Calculating effect of current displacement on resistance
                mu_0 = 4*pi*1e-7;
                
                zeta_Rs = sqrt((pi * mu_0 * f_elec) ./ obj.rho_s) .* obj.h_slot;
                K_Rs = zeta_Rs .* (sinh(2*zeta_Rs) + sin(2*zeta_Rs)) ./ (cosh(2*zeta_Rs) - cos(2*zeta_Rs));
                
                zeta_Rr = sqrt((pi*mu_0*obj.f)/obj.rho_s) * D_rotorbar; % Assumes circular cross section of rotor bars
                K_Rr = zeta_Rr .* (sinh(2*zeta_Rr) + sin(2*zeta_Rr)) ./ (cosh(2*zeta_Rr) - cos(2*zeta_Rr));
            end
            Rs_AC = obj.Rs_DC .* K_Rs;
            Rr_AC = obj.Rr_DC*K_Rr;


            % 2. Wrap angle to 0:2pi
            theta_wrapped = mod(Theta, 2*pi);
            
            % 3. Calculate LUT index as described in Lipo paper (Linear interpolation)
            % exact position in the array (1-based)
            pos_float = (theta_wrapped / obj.d_theta) + 1; 
            
            idx_1 = floor(pos_float);
            idx_2 = idx_1 + 1;
            
            % Handle Wrap-around at 360 deg
            if idx_1 >= obj.Num_Steps
                idx_1 = obj.Num_Steps;
                idx_2 = 1; % Wrap to start
            end
            
            % Interpolation weight
            w = pos_float - idx_1; 
            
            % 4. Fetch matricies
            % Interpolate Inverse L
            InvL = (1-w) * obj.InvL_LUT(:,:,idx_1) + w * obj.InvL_LUT(:,:,idx_2);
            
            % Interpolate G (dL/dTheta)
            G = (1-w) * obj.G_LUT(:,:,idx_1)+ w * obj.G_LUT(:,:,idx_2);

            % Dynamic Resistance Matrix
            R_vec = [repmat(Rs_AC, Qs, 1); repmat(Rr_AC, Qr, 1)];
            obj.R_Matrix = diag(R_vec);
            
            % 5. Solve Dynamics
            I = obj.Currents;
            V = [V_stator; zeros(Qr, 1)];
            
            % dI/dt = InvL * (V - R*I - w*G*I)
            % Note: G*I is the back-EMF flux derivative term
            BackEMF_Voltage = Omega * (G * I);
            Resistive_Voltage = obj.R_Matrix * I;
            
            Voltage_Sum = V - Resistive_Voltage - BackEMF_Voltage;
            
            dIdt = InvL * Voltage_Sum;
            
            % 6. Update State
            I_new = I + dIdt * obj.Ts;
            obj.Currents = I_new;
            
            % 7. Output results
            Torque = 0.5 * I' * G * I;
            I_stator = I_new(1:Qs);
        end

        function resetImpl(obj)
            % Reset with sizes from mparams if available, else property defaults
            % Easier to just query current state size
            obj.Currents = zeros(size(obj.Currents));
        end
        
        %% Discrete State Specification
        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, name)
            % Specify size, data type, and complexity of discrete state
            if strcmp(name, 'Currents')
                sz = [obj.Qs + obj.Qr, 1];  % [110 x 1]
                dt = 'double';
                cp = false;  % real, not complex
            end
        end

        %% Define output sizes and types
        function [sz1, sz2] = getOutputSizeImpl(obj)
            % Output 1: I_stator is [Qs x 1] = [60 x 1]
            sz1 = [obj.Qs, 1];
            % Output 2: Torque is scalar [1 x 1]
            sz2 = [1, 1];
        end

        function [dt1, dt2] = getOutputDataTypeImpl(~) 
            dt1 = 'double';
            dt2 = 'double';
        end
        
        function [cp1, cp2] = isOutputComplexImpl(~)
            cp1 = false;
            cp2 = false;
        end

        function [fz1, fz2] = isOutputFixedSizeImpl(~)
            % Outputs are fixed size (not variable)
            fz1 = true;
            fz2 = true;
        end
    end
end

% Helper function for GCD of floating point numbers
function result = gcd_float(a, b)
    while b > 1e-5 % Tolerance
        temp = b;
        b = mod(a, b);
        a = temp;
    end
    result = a;
end