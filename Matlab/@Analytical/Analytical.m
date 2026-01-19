classdef Analytical < dynamicprops
    methods (Static)
        function obj = LoadFromJSON(filename)
            % fromFile  Create Analytical by loading JSON file
            s = jsondecode(fileread(filename));
            obj = Analytical();
            obj.addPropsFromStruct(s);
        end
    end

    methods
        function addPropsFromStruct(obj, s)
            % addPropsFromStruct  Recursively add properties from struct s
            if ~isstruct(s)
                error('Input must be a struct (result of jsondecode).');
            end
            fields = fieldnames(s);
            for k = 1:numel(fields)
                name = fields{k};
                val  = s.(name);

                % Ensure valid property name
                if ~isvarname(name)
                    error('Field name "%s" is not a valid MATLAB identifier.', name);
                end

                % Add dynamic property
                if isempty(findprop(obj, name))
                    addprop(obj, name);
                end

                % If value is a struct, create nested Analytical
                if isstruct(val)
                    child = Analytical();
                    child.addPropsFromStruct(val);
                    obj.(name) = child;
                else
                    % All values assumed numeric/double in JSON: ensure double
                    if isnumeric(val)
                        obj.(name) = double(val);
                    else
                        % Optionally attempt conversion for scalar-like values
                        num = str2double(val);
                        if ~isnan(num)
                            obj.(name) = double(num);
                        else
                            % If non-numeric encountered, store as-is
                            obj.(name) = val;
                        end
                    end
                end
            end
        end
    end
end
