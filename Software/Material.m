% Stores material data read from material_data.txt for a material config.
% Initialized with a configName input, which it uses to search the config
% file for what data to pull.
% -------------------------------------------------------------------------
% Dependencies
%   1) material_data.txt
%   2) GetConfigData.m Version 2.1
% -------------------------------------------------------------------------
% Assumptions
%   1) All units are SI
% -------------------------------------------------------------------------
% Comments
%   1) Vectorized. Pass in configName as an array and a same-sized array of
%       Material objects will be returned.
%   2) getK and getSigmaY are vectorized over T of any size.
% -------------------------------------------------------------------------
% Nomenclature
%   T = temperature (K)
% -------------------------------------------------------------------------
% Document Version 1.0, earlier versions:
%   - <Version>
% -------------------------------------------------------------------------
% MATLAB Version R2025b, also compatible with:
%   - R2024a
% -------------------------------------------------------------------------
% Developed by Alex Vance (AlexVance00 on Github)
classdef Material

    properties (Constant, Hidden)
        configFile = "material_data.txt";
    end

    properties (SetAccess = private)
        name % Config name
        rho % Density (kg/m^3)
        cp % Specific heat at constant pressure (J/kg-K)
        T_melt % Melting temperature (K)
        E % Elastic modulus (Pa)
        CTE % Coefficient of thermal expansion (1/K)
        epsilon % Surface roughness (m)
        k % thermal conductivity vs. temperature, W/(m-K) vs. K
        sigma_y % Yield stress vs. temperature, Pa vs. K
    end

    methods (Access = public)

        % Constructor
        function obj = Material(configName)
            % If no arguments were passed, just return empty object
            if nargin == 0
                return;
            end

            % Check if configName was correctly passed as string
            if ~isstring(configName)
                % Check if it was a char instead
                if ~ischar(configName)
                    error("Passed ""configName"" argument must be string or char types, but was passed as %s", class(configName));
                else
                    configName = string(configName);
                end
            end

            data = GetConfigData(configFile = Material.configFile, configName = configName);

            % If numConfigNames is not 1, flag there as being multiple and
            % preallocate array of obj types for speed's sake
            numConfigNames = numel(configName);
            flagMultipleConfigNames = false;
            if numConfigNames ~= 1
                flagMultipleConfigNames = true;
                obj(numConfigNames) = Material();
            end

            % Assign data to properties
            for i = 1:numConfigNames
                thisData = data(i);
                theseFieldnames = string(fieldnames(thisData));
                for j = 1:numel(theseFieldnames)
                    obj(i).(theseFieldnames(j)) = thisData.(theseFieldnames(j));
                end
                obj(i).name = configName(i);
            end

            % Make sure array of obj type is same size as configName, if
            % it's not 1
            if flagMultipleConfigNames
                obj = reshape(obj, size(configName));
            end

            return;
        end

        function k = getK(obj, T)
        % Gets thermal conductivity at temperature T in W(/m-K)
        % -----------------------------------------------------------------
        % Arguments
        %   T = temperature (K), any size
        % -----------------------------------------------------------------
        % Assumptions
        %   1) For T values outside the table's range, an error is thrown
        % -----------------------------------------------------------------
            k = obj.evalProperty("k", T);
        end

        function sigma_y = getSigmaY(obj, T)
        % Gets yield stress at temperature T in Pa
        % -----------------------------------------------------------------
        % Arguments
        %   T = temperature (K), any size
        % -----------------------------------------------------------------
        % Assumptions
        %   1) For T values outside the table's range, an error is thrown
        % -----------------------------------------------------------------
            sigma_y = obj.evalProperty("sigma_y", T);
        end

    end

    methods (Static, Access = public)

        function printValidMaterials()
        % Prints the names of all configs. in material_data.txt
            Material.listConfigNames();
        end

    end

    methods (Access = private)

        function value = evalProperty(obj, propName, T)
        % Evaluates a property that is either a constant or a correlation,
        % returning an array the same size as T
            prop = obj.(propName);
            if isempty(prop)
                error("Material ""%s"" has no ""%s"" in %s", obj.name, propName, Material.configFile);
            end
            value = interp1(prop(:, 1), prop(:, 2), T);
        end

    end

    methods (Static, Access = private)

        function data = collapseCorrelations(data)
        % Turns GetConfigData's flags into one struct(model, coeffs) per
        % correlation property, and removes the <property>_model and
        % <property>_coeffs fields.
            fieldNames = string(fieldnames(data));
            for i = 1:numel(fieldNames)
                fieldName = fieldNames(i);
                % Earlier correlations in this loop remove their own _model and
                % _coeffs fields, which are still in fieldNames
                if ~isfield(data, fieldName)
                    continue;
                end
                value = data.(fieldName);
                if isstring(value) && isscalar(value) && value == "flagCorrelation"
                    modelName = fieldName + "_model";
                    coeffsName = fieldName + "_coeffs";
                    if ~isfield(data, modelName) || ~isfield(data, coeffsName)
                        error("Property ""%s"" is ""correlation"" in %s but ""%s"" and/or ""%s"" are missing", fieldName, Material.configFile, modelName, coeffsName);
                    end
                    data.(fieldName) = struct("model", data.(modelName), "coeffs", data.(coeffsName));
                    data = rmfield(data, [modelName, coeffsName]);
                end
            end
        end

        function listConfigNames()
        % Prints config names found in the config file
            configFileSearchMatches = dir(fullfile(pwd, "**", Material.configFile));
            if isempty(configFileSearchMatches)
                error("Config file ""%s"" not found in working directory ""%s""\n", Material.configFile, pwd);
            end
            rawText = readlines(fullfile(configFileSearchMatches(1).folder, configFileSearchMatches(1).name));
            allConfigsText = rawText(find(rawText == "%", 1) + 1:end);
            isHeaderLine = startsWith(allConfigsText, "\") & endsWith(allConfigsText, "\") & allConfigsText ~= "\END\";
            names = erase(allConfigsText(isHeaderLine), "\");
            fprintf("Valid material names in %s:\n", Material.configFile);
            fprintf("  %s\n", names);
        end

    end

end
