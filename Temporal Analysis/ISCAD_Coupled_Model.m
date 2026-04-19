classdef ISCAD_Coupled_Model < matlab.System
    % ISCAD_Coupled_Model
    % Coupled Circuit Model (Phase Variable Model) for ISCAD machine.
    % Uses pre-calculated LUTs for L(theta) and G(theta) = dL/dtheta.
    % Grid points aligned to stator/rotor slotting via GCD of slot pitches.
    %
    % CONVENTIONS:
    %   Theta  = mechanical rotor angle [rad] (0 to 2*pi per revolution)
    %   Omega  = mechanical angular velocity [rad/s]
    %   p      = number of poles (not pole pairs)
    %   V_stator = [Qs x 1] stator voltage vector
    %   I_feedback = [Qs+Qr x 1] full current state (stator + rotor)

    % Public Properties (set before simulation, cannot change during run)
    properties (Nontunable)
        mparams = struct();
        % LUT Resolution Multiplier
        % 1 = Points exactly at slot alignments. 10 = 10 points between alignments.
        Resolution_Factor = 10;
        f_sw = 20e3; % Switching Frequency [Hz]
        Qs = 60; % Number of Stator Slots
        Qr = 70; % Number of Rotor Bars (from JSON)
        p = 2; % Number of Poles (not pole pairs)

        rho_s = 0; % Resistivity of Conductor [Ohm*m]
        h_slot = 0; % Slot Height [m]
        ag = 0; % Air Gap [m]
        l_stack = 0; % Stack Length [m]
        w_slot = 0; % Slot Width [m] (may be [inner, outer])
        r = 0; % Bore Radius [m]
        h_so = 0; % Slot Opening Height [m]
        w_so = 0; % Slot Opening Width [m]

        Ts = 1e-6;  % Step time [s]
    end

    properties (SetAccess = protected, GetAccess = public)
        % Derived parameters
        Rs_DC       % Stator DC resistance [Ohm]
        Rr_DC       % Rotor DC resistance [Ohm]
        D_rotorbar  % Equivalent rotor bar diameter [m]
        L_ss        % Per-slot magnetising self-inductance [H]
        L_ls        % Stator leakage inductance (DC) [H]
        L_lr        % Rotor leakage inductance (DC) [H]

        L_LUT       % 3D Matrix: [N_sys x N_sys x Steps] Stores L(theta)
        G_LUT       % 3D Matrix: [N_sys x N_sys x Steps] Stores dL/dTheta
        R_Matrix    % Resistance matrix (updated at runtime)

        d_theta     % Step size of the LUT [rad]
        Num_Steps   % Number of angular steps in 360 degrees
    end

    methods(Access = protected)

        function setupImpl(obj)
            % setupImpl: Pre-calculate inductance LUTs over one mechanical
            % revolution. Called once before simulation starts.

            mu_0 = 4*pi*1e-7;

            % --- Calculate DC Resistances ---
            % Stator (w_slot may be [inner, outer]; use index 1)
            obj.Rs_DC = obj.rho_s * obj.l_stack / (obj.h_slot * obj.w_slot(1));

            % Rotor (total rotor copper area scaled from stator)
            A_rotorbar = obj.h_slot * obj.w_slot(1) * obj.Qs / obj.Qr;
            obj.Rr_DC = obj.rho_s * obj.l_stack / A_rotorbar;
            obj.D_rotorbar = sqrt(4*A_rotorbar/pi); % Equivalent diameter

            % --- Calculate Inductances ---
            % Per-slot magnetising self-inductance (smooth-bore approximation)
            obj.L_ss = (pi/2) * mu_0 * (obj.r * obj.l_stack / obj.ag);

            % Leakage inductance at DC (no skin-effect reduction).
            % Skin effect is already captured in the dynamic resistance R(f).
            % Using DC leakage for the inductance matrix improves numerical
            % conditioning (~6,700 vs ~77,000 with switching-frequency leakage).
            % TODO: Update h_so and w_so with actual slot opening dimensions
            lambda_slot = (obj.h_slot / (3*obj.w_slot(1))) + (obj.h_so / obj.w_so);
            obj.L_ls = mu_0 * obj.l_stack * lambda_slot;
            obj.L_lr = obj.L_ls; % Assume equal leakage (improve when rotor geometry is defined)

            % Local copies for LUT generation loop
            L_ss_val = obj.L_ss;
            L_ls_val = obj.L_ls;
            L_lr_val = obj.L_lr;

            fprintf('L_ss (magnetising) = %.3e H, L_ls (leakage DC) = %.3e H\n', ...
                    L_ss_val, L_ls_val);

            % --- Determine Optimal LUT Discretisation ---
            % Slot pitches in degrees
            deg_s = 360 / obj.Qs;
            deg_r = 360 / obj.Qr;

            % GCD ensures grid hits stator-rotor alignment points
            common_pitch = gcd_float(deg_s, deg_r);
            step_deg = common_pitch / obj.Resolution_Factor;
            obj.Num_Steps = round(360 / step_deg);
            obj.d_theta = deg2rad(step_deg);

            fprintf('LUT: %g angular steps (%.4f deg each)\n', ...
                    obj.Num_Steps, step_deg);

            % --- Pre-Calculate LUTs ---
            N_sys = obj.Qs + obj.Qr;
            obj.L_LUT = zeros(N_sys, N_sys, obj.Num_Steps);
            obj.G_LUT = zeros(N_sys, N_sys, obj.Num_Steps);

            % Slot angle vectors (mechanical angles, 0-indexed)
            as_base = (0:obj.Qs-1)' * (2*pi/obj.Qs);
            ar_base = (0:obj.Qr-1)' * (2*pi/obj.Qr);

            worst_cond = 0;

            for k = 1:obj.Num_Steps
                theta = (k-1) * obj.d_theta;
                ar_grid = ar_base + theta;

                % -- Stator-Stator block (independent of theta) --
                % Triangular mutual inductance from integration of square MMF
                diff_ss = mod(as_base - as_base' + pi, 2*pi) - pi;
                Lss = L_ss_val * (1 - 2*abs(diff_ss)/pi);
                % Add stator leakage to diagonal
                Lss(1:obj.Qs+1:end) = Lss(1:obj.Qs+1:end) + L_ls_val;

                % -- Rotor-Rotor block (independent of theta; theta cancels) --
                diff_rr = mod(ar_grid - ar_grid' + pi, 2*pi) - pi;
                Lrr = L_ss_val * (1 - 2*abs(diff_rr)/pi);
                % Add rotor leakage to diagonal
                Lrr(1:obj.Qr+1:end) = Lrr(1:obj.Qr+1:end) + L_lr_val;

                % -- Stator-Rotor mutual (varies with theta) --
                diff_sr = mod(as_base - ar_grid' + pi, 2*pi) - pi;
                Msr = L_ss_val * (1 - 2*abs(diff_sr)/pi);

                % Assemble full inductance matrix
                L_mat = [Lss, Msr; Msr', Lrr];

                % -- dL/dtheta: derivative of Msr w.r.t. mechanical angle --
                % d/dtheta of triangle wave = square wave
                % Chain rule: d|as-ar-theta|/dtheta = -sign(as-ar-theta)
                % So dMsr/dtheta = L_ss * (2/pi) * sign(diff_sr)
                dMsr = L_ss_val * (2/pi) * sign(diff_sr);
                % Lss and Lrr blocks are theta-independent, so their
                % derivatives are zero
                G_mat = [zeros(obj.Qs), dMsr; dMsr', zeros(obj.Qr)];

                % Store L (not inv(L)) to avoid interpolation-of-inverse error
                obj.L_LUT(:,:,k) = L_mat;
                obj.G_LUT(:,:,k) = G_mat;

                % Track worst condition number for diagnostics
                c = cond(L_mat);
                if c > worst_cond
                    worst_cond = c;
                end
            end

            fprintf('LUT pre-calculation done. Worst condition number: %.0f\n', worst_cond);
            if worst_cond > 10000
                warning('ISCAD:ConditionNumber', ...
                    'L matrix condition number %.0f is high. Check leakage inductance values.', ...
                    worst_cond);
            end
        end

        function [dIdt, Torque] = stepImpl(obj, V_stator, I_feedback, Theta, Omega, f_elec)
            % stepImpl: Called every simulation timestep.
            %   V_stator  [Qs x 1]     - Stator voltage vector (from PWM)
            %   I_feedback [Qs+Qr x 1] - Full current state (from integrator)
            %   Theta     [scalar]      - Mechanical rotor angle [rad]
            %   Omega     [scalar]      - Mechanical angular velocity [rad/s]
            %   f_elec    [scalar]      - Electrical frequency [Hz] (for skin effect)

            Qs = obj.Qs;
            Qr = obj.Qr;

            % 1. Dynamic AC resistance (skin effect at fundamental frequency)
            if f_elec < 1
                K_Rs = 1;
                K_Rr = 1;
            else
                mu_0 = 4*pi*1e-7;

                zeta_Rs = sqrt((pi * mu_0 * f_elec) / obj.rho_s) * obj.h_slot;
                K_Rs = zeta_Rs * (sinh(2*zeta_Rs) + sin(2*zeta_Rs)) / (cosh(2*zeta_Rs) - cos(2*zeta_Rs));

                zeta_Rr = sqrt((pi * mu_0 * f_elec) / obj.rho_s) * obj.D_rotorbar;
                K_Rr = zeta_Rr * (sinh(2*zeta_Rr) + sin(2*zeta_Rr)) / (cosh(2*zeta_Rr) - cos(2*zeta_Rr));
            end
            Rs_AC = obj.Rs_DC * K_Rs;
            Rr_AC = obj.Rr_DC * K_Rr;

            % 2. Wrap mechanical angle to [0, 2*pi)
            theta_wrapped = mod(Theta, 2*pi);

            % 3. LUT interpolation (linear, as per Lipo methodology)
            pos_float = (theta_wrapped / obj.d_theta) + 1; % 1-based index
            idx_1 = floor(pos_float);
            idx_2 = idx_1 + 1;

            % Handle wrap-around at 360 deg
            if idx_1 >= obj.Num_Steps
                idx_1 = obj.Num_Steps;
                idx_2 = 1;
            end

            % Interpolation weight
            w = pos_float - idx_1;

            % 4. Interpolate L and G matrices
            L_interp = (1-w) * obj.L_LUT(:,:,idx_1) + w * obj.L_LUT(:,:,idx_2);
            G = (1-w) * obj.G_LUT(:,:,idx_1) + w * obj.G_LUT(:,:,idx_2);

            % 5. Build resistance matrix
            R_vec = [repmat(Rs_AC, Qs, 1); repmat(Rr_AC, Qr, 1)];

            % 6. Solve state-space ODE
            % V = R*I + L*dI/dt + Omega*G*I
            % => dI/dt = L \ (V - R*I - Omega*G*I)
            I = I_feedback;
            V = [V_stator; zeros(Qr, 1)]; % Rotor is short-circuited

            BackEMF_Voltage = Omega * (G * I);
            Resistive_Voltage = R_vec .* I; % Element-wise (avoids building full diag matrix)

            rhs = V - Resistive_Voltage - BackEMF_Voltage;
            dIdt = L_interp \ rhs; % Backslash solve (numerically stable)

            % 7. Electromagnetic torque
            % T_e = (p/2) * (1/2) * I' * dL/dtheta_mech * I
            Torque = (obj.p / 2) * 0.5 * (I' * G * I);
        end


        %% Define output sizes and types
        function [sz1, sz2] = getOutputSizeImpl(obj)
            sz1 = [obj.Qs + obj.Qr, 1]; % dIdt
            sz2 = [1, 1];               % Torque (scalar)
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
            fz1 = true;
            fz2 = true;
        end
    end
end

% Helper function for GCD of floating point numbers
function result = gcd_float(a, b)
    while b > 1e-5
        temp = b;
        b = mod(a, b);
        a = temp;
    end
    result = a;
end
