classdef ISCAD_Coupled_Model_3ph < matlab.System
    % ISCAD_Coupled_Model_3ph
    % 3-Phase Coupled Circuit Model for ISCAD test rig.
    % Models 3 conductor bars on a shared magnetic core with mutual
    % coupling via triangular mutual inductance waveform.
    % No rotor — purely stator-side harmonic coupling investigation.
    %
    % This model is designed to match a physical 3-phase test rig where
    % 3 bars share a common magnetic core and are each driven by an
    % independent half-bridge sharing a DC bus.
    %
    % CONVENTIONS:
    %   V_stator = [3 x 1] voltage vector (from 3 half-bridges)
    %   I_feedback = [3 x 1] current state (from integrator)
    %
    % The governing equation (no rotor, no rotation):
    %   V = R*I + L*dI/dt
    %   dI/dt = L \ (V - R*I)

    properties (Nontunable)
        % Self inductance of each bar [H]
        L_self = 0.25e-3;
        % Mutual inductance between adjacent bars [H]
        % For 3 bars at 120 deg spacing on ISCAD geometry:
        % M = L_ss * (1 - 2*|120deg|/pi) = L_ss * (1 - 4/3) = -L_ss/3
        % Set to actual measured/calculated value for test rig
        L_mutual_12 = 0;
        L_mutual_13 = 0;
        L_mutual_23 = 0;
        % Phase resistance [Ohm]
        R_ph = 0.8e-3;
        % If true, all mutual inductances are set equal to L_mutual_12
        Use_symmetric_mutual = true;
    end

    properties (SetAccess = protected, GetAccess = public)
        L_mat   % 3x3 inductance matrix
        R_vec   % 3x1 resistance vector
    end

    methods(Access = protected)

        function setupImpl(obj)
            % Build the 3x3 inductance matrix
            if obj.Use_symmetric_mutual
                M = obj.L_mutual_12;
                obj.L_mat = [obj.L_self, M,              M;
                             M,          obj.L_self,      M;
                             M,          M,               obj.L_self];
            else
                obj.L_mat = [obj.L_self,      obj.L_mutual_12, obj.L_mutual_13;
                             obj.L_mutual_12, obj.L_self,      obj.L_mutual_23;
                             obj.L_mutual_13, obj.L_mutual_23, obj.L_self];
            end

            obj.R_vec = [obj.R_ph; obj.R_ph; obj.R_ph];

            fprintf('\n--- 3-Phase Test Rig Model ---\n');
            fprintf('L_self = %.3e H, R_ph = %.3e Ohm\n', obj.L_self, obj.R_ph);
            fprintf('L matrix condition number: %.1f\n', cond(obj.L_mat));

            if any(eig(obj.L_mat) <= 0)
                error('ISCAD_3ph:NotPositiveDefinite', ...
                    'Inductance matrix is not positive definite. Check mutual inductance values.');
            end
        end

        function [dIdt, I_sum] = stepImpl(obj, V_stator, I_feedback)
            % stepImpl: Called every simulation timestep.
            %   V_stator   [3 x 1] - Phase voltages from half-bridges
            %   I_feedback [3 x 1] - Current state from integrator
            %
            % Outputs:
            %   dIdt  [3 x 1] - Current derivative for integrator
            %   I_sum [1 x 1] - Sum of all phase currents (should be ~0
            %                    for balanced operation; nonzero indicates
            %                    common-mode / harmonic coupling)

            I = I_feedback;

            % V = R*I + L*dI/dt  =>  dI/dt = L \ (V - R*I)
            Resistive_Voltage = obj.R_vec .* I;
            rhs = V_stator - Resistive_Voltage;
            dIdt = obj.L_mat \ rhs;

            % Common-mode current diagnostic
            I_sum = I(1) + I(2) + I(3);
        end

        %% Define output sizes and types
        function [sz1, sz2] = getOutputSizeImpl(~)
            sz1 = [3, 1]; % dIdt
            sz2 = [1, 1]; % I_sum (scalar)
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
