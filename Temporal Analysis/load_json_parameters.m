function mparams = load_json_parameters(json_file)
% load_json_parameters - Load parameters from JSON file into flattened struct
%
% Usage:
%   mparams = load_json_parameters('ISCAD_parameters.json')
%
% This function reads the JSON parameter file and flattens the hierarchy,
% creating a single structure 'mparams' containing all parameters.
% The struct is returned and also assigned to the 'base' workspace
% for Simulink access.

    % Read and parse JSON file
    fid = fopen(json_file, 'r');
    if fid == -1
        error('Cannot open file: %s', json_file);
    end
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    
    raw_params = jsondecode(str);
    fprintf('Parameters loaded from: %s\n', json_file);
    
    % Initialize flattened struct
    mparams = struct();
    
    % Recursively flatten the structure
    fields = fieldnames(raw_params);
    for i = 1:length(fields)
        val = raw_params.(fields{i});
        if isstruct(val)
            subfields = fieldnames(val);
            for j = 1:length(subfields)
                % Check for name collision (e.g. Apk in multiple places)
                % We prioritize the first occurrence or specific logic if needed
                % For now, simpler override logic (last one wins if collision)
                subval = val.(subfields{j});
                
                % Special handling if needed (e.g. rename)
                mparams.(subfields{j}) = subval;
            end
        else
            mparams.(fields{i}) = val;
        end
    end
    
    % Calculate derived parameters
    mparams.m = mparams.Qs / mparams.p;
    mparams.r = mparams.r; % Ensure r is available by name
    mparams.r_r = mparams.r - mparams.ag;
    
    % Assign to base workspace for Simulink
    assignin('base', 'mparams', mparams);
end
