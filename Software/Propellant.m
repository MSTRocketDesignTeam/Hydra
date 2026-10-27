% Stores propellant data read from propellant_data.txt for a propellant
% config. Initialized with a configName input, which it uses to search the
% config file for what data to pull. Replaces the old abstract Propellant
% superclass + one-subclass-per-propellant scheme: every propellant is now a
% Propellant object, and any property that is not a constant (density,
% conductivity, viscosity, vapor pressure, ...) is evaluated by a correlation
% *model* selected in the config file, with its coefficients also stored in
% the config file. Models are written once, here, and shared by every
% propellant that uses them.
% -------------------------------------------------------------------------
% Dependencies
%   1) propellant_data.txt
%   2) GetConfigData.m Version 2.1
% -------------------------------------------------------------------------
% Assumptions
%   1) A property in the config is either a constant number, or
%       "correlation". A "correlation" property <p> must have <p>_model
%       and <p>_coeffs lines in the config. Instance data size stays small:
%       a correlation property is stored as struct(model, coeffs), where
%       coeffs is a short numeric array, not a function handle or table.
%   2) Tabulated (.xlsx) properties are not supported by the evaluation
%       methods (the "table" model, with the table stored in the config
%       file itself, replaces them).
%   3) All units are SI.
% -------------------------------------------------------------------------
% Comments
%   1) Vectorized. Pass in configName as an array and a same-sized array of
%       Propellant objects will be returned.
%   2) Property evaluation methods (getRho, getCp, getK, getMu, getPvap) are
%       vectorized over T (and rho) of any, mutually compatible, size.
%   3) To add a new correlation model, add a case to evalModel and a private
%       static method for it, and document its coefficient layout in both
%       that method's header and propellant_data.txt's instructions.
% -------------------------------------------------------------------------
% Nomenclature
%   T = temperature (K)
%   rho = density (kg/m^3)
%   Tr = reduced temperature, T/T_crit
%   rho_r = reduced density, rho/rho_crit
% -------------------------------------------------------------------------
% Document Version 2.0, earlier versions:
%   - 1.0 (abstract superclass with Ethanol and N2O subclasses)
% -------------------------------------------------------------------------
% MATLAB Version R2024b, also compatible with:
%   - R2025b
% -------------------------------------------------------------------------
% Developed by Alex Vance (AlexVance00 on Github)
classdef Propellant

    properties (Constant, Hidden)
        NA = 6.02214076E23; % Avogadro's Number
        R_universal = 8.314; % Universal gas constant
        configFile = "propellant_data.txt";
    end

    properties (SetAccess = private)
        name % Config name, also the CEA propellant name
        phase % "LIQUID" or "GAS"
        type % "FUEL" or "OXIDIZER"
        cost % Cost per unit mass (USD/kg)
        mlr_wgt % Molar weight (kg/mol)
        T_crit % Critical temperature (K)
        P_crit % Critical pressure (Pa)
        rho_crit % Critical density (kg/m^3)
        % The following are either a constant number or a struct with
        % fields model and coeffs (see Assumptions). Use the get* methods
        % to evaluate them, they work either way.
        rho % Density (kg/m^3)
        cp % Specific heat at constant pressure (J/kg-K)
        k % Thermal conductivity (W/m-K)
        mu % Dynamic viscosity (Pa-s)
        p_vap % Vapor pressure (Pa)
    end

    methods (Access = public)

        % Constructor
        function obj = Propellant(configName)
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

            data = GetConfigData(configFile = Propellant.configFile, configName = configName);

            % If numConfigNames is not 1, flag there as being multiple and
            % preallocate array of obj types for speed's sake
            numConfigNames = numel(configName);
            flagMultipleConfigNames = false;
            if numConfigNames ~= 1
                flagMultipleConfigNames = true;
                obj(numConfigNames) = Propellant();
            end

            % Not every config has every property, so only assign what was
            % in the config. Anything missing stays empty and the get*
            % methods will say so if it's asked for.
            for i = 1:numConfigNames
                thisData = Propellant.collapseCorrelations(data(i));
                dataFieldNames = string(fieldnames(thisData));
                for j = 1:numel(dataFieldNames)
                    obj(i).(dataFieldNames(j)) = thisData.(dataFieldNames(j));
                end
            end

            % Make sure array of obj type is same size as configName, if
            % it's not 1
            if flagMultipleConfigNames
                obj = reshape(obj, size(configName));
            end

            return;
        end

        function rho = getRho(obj, T)
        % Density at temperature T                          rho in (kg/m^3)
        % -----------------------------------------------------------------
        % Arguments
        %   T = temperature (K), any size
        % -----------------------------------------------------------------
            rho = obj.evalProperty("rho", T, []);
        end

        function cp = getCp(obj, T)
        % Specific heat at constant pressure at temperature T
        %                                                   cp in (J/kg-K)
        % -----------------------------------------------------------------
        % Arguments
        %   T = temperature (K), any size
        % -----------------------------------------------------------------
            cp = obj.evalProperty("cp", T, []);
        end

        function k = getK(obj, T, rho)
        % Thermal conductivity at temperature T and density rho
        %                                                    k in (W/m-K)
        % -----------------------------------------------------------------
        % Arguments
        %   T = temperature (K), any size
        %   rho = density (kg/m^3), same size as T (or scalar). Only needed
        %       if the config's model for k depends on density
        % -----------------------------------------------------------------
            if nargin < 3
                rho = [];
            end
            k = obj.evalProperty("k", T, rho);
        end

        function mu = getMu(obj, T, rho)
        % Dynamic viscosity at temperature T and density rho   mu in (Pa-s)
        % -----------------------------------------------------------------
        % Arguments
        %   T = temperature (K), any size
        %   rho = density (kg/m^3), same size as T (or scalar). Only needed
        %       if the config's model for mu depends on density
        % -----------------------------------------------------------------
            if nargin < 3
                rho = [];
            end
            mu = obj.evalProperty("mu", T, rho);
        end

        function p_vap = getPvap(obj, T)
        % Vapor pressure at temperature T                     p_vap in (Pa)
        % -----------------------------------------------------------------
        % Arguments
        %   T = temperature (K), any size
        % -----------------------------------------------------------------
            p_vap = obj.evalProperty("p_vap", T, []);
        end

    end

    methods (Static, Access = public)

        function printValidNames()
        % Prints the names of all configs. in propellant_data.txt
            Propellant.listConfigNames();
        end

    end

    methods (Access = private)

        function value = evalProperty(obj, propName, T, rho)
        % Evaluates a property that is either a constant or a correlation,
        % returning an array the same size as T
            prop = obj.(propName);
            if isempty(prop)
                error("Propellant ""%s"" has no ""%s"" in %s", obj.name, propName, Propellant.configFile);
            end
            if isnumeric(prop)
                % Constant property
                value = prop .* ones(size(T));
            elseif isstruct(prop) && isfield(prop, "model")
                value = obj.evalModel(prop.model, prop.coeffs, T, rho);
            else
                error("Property ""%s"" of ""%s"" is neither a constant nor a correlation; tabulated (.xlsx) properties aren't supported by the get* methods", propName, obj.name);
            end
        end

        function value = evalModel(obj, model, coeffs, T, rho)
        % Dispatches to the correlation model named in the config
            switch string(model)
                case "rackett"
                    value = obj.modelRackett(coeffs, T);
                case "assael_ethanol"
                    if isempty(rho)
                        error("Model ""assael_ethanol"" needs a density input");
                    end
                    value = obj.modelAssaelEthanolK(coeffs, T, rho);
                case "sotiriadou_ethanol"
                    if isempty(rho)
                        error("Model ""sotiriadou_ethanol"" needs a density input");
                    end
                    value = obj.modelSotiriadouEthanolMu(coeffs, T, rho);
                case "antoine"
                    value = Propellant.modelAntoine(coeffs, T);
                case "table"
                    value = interp1(coeffs(:, 1), coeffs(:, 2), T);
                otherwise
                    error("Unknown correlation model ""%s"" for propellant ""%s""", model, obj.name);
            end
        end

        function rho = modelRackett(obj, coeffs, T)
        % Saturated liquid density by the Rackett equation   rho in (kg/m^3)
        % -----------------------------------------------------------------
        % Arguments
        %   coeffs = [Zc]
        %       Zc = compressibility factor used by the correlation (~)
        %   T = temperature (K)
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Critical volume is derived from the ideal-gas-like relation
        %       Vc = Zc * R * T_crit / P_crit, using this object's T_crit
        %       and P_crit, and the same Zc is used for the exponent base.
        % -----------------------------------------------------------------
        % Sources
        %   1) Rackett, H.G. (1970) "Equation of State for Saturated
        %       Liquids", J. Chem. Eng. Data, 15(4), 514-517.
        % -----------------------------------------------------------------
            Zc = coeffs(1);
            Vc = Zc .* Propellant.R_universal .* obj.T_crit ./ obj.P_crit;
            Tr = T ./ obj.T_crit;
            rho = obj.mlr_wgt ./ (Vc .* Zc.^((1 - Tr).^(2 ./ 7)));
        end

        function k = modelAssaelEthanolK(obj, coeffs, T, rho)
        % Thermal conductivity of ethanol, Assael et al. 2010   k in (W/m-K)
        % -----------------------------------------------------------------
        % Arguments
        %   coeffs = 22 element row vector:
        %       (1:6)   lambda0 numerator, coefficients of Tr^0 ... Tr^5
        %       (7:9)   lambda0 denominator, coefficients of Tr^0 ... Tr^2
        %       (10:14) B1_i, i = 1 ... 5
        %       (15:19) B2_i, i = 1 ... 5
        %       (20:22) critical enhancement [amplitude, offset, exponent
        %               coefficient]
        %   T = temperature (K)
        %   rho = density (kg/m^3)
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Dilute-gas term is in mW/m-K in the paper, converted to
        %       W/m-K here.
        % -----------------------------------------------------------------
        % Sources
        %   1) Assael, M.J., et al., Reference Data for the Thermal
        %       Conductivity of Ethanol, J. Phys. Chem. Ref. Data, 2010
        % -----------------------------------------------------------------
            Tr = T ./ obj.T_crit;
            rho_r = rho ./ obj.rho_crit;

            % Dilute-gas term, ratio of polynomials in Tr (summed term by
            % term so it reads like the paper and works for any size of T)
            numerator = zeros(size(Tr));
            for i = 1:6
                numerator = numerator + coeffs(i) .* Tr.^(i - 1);
            end
            denominator = zeros(size(Tr));
            for i = 1:3
                denominator = denominator + coeffs(6 + i) .* Tr.^(i - 1);
            end
            lambda0 = numerator ./ denominator .* 1E-3;

            % Residual term, sum over i = 1 ... 5
            B1 = coeffs(10:14);
            B2 = coeffs(15:19);
            delta_lambda = zeros(size(Tr .* rho_r));
            for i = 1:5
                delta_lambda = delta_lambda + (B1(i) + B2(i) .* Tr) .* rho_r.^i;
            end

            % Critical enhancement
            delta_lambda_c = coeffs(20) ./ (coeffs(21) + abs(Tr - 1)) .* exp(-(coeffs(22) .* (rho_r - 1)).^2);

            k = lambda0 + delta_lambda + delta_lambda_c;
        end

        function mu = modelSotiriadouEthanolMu(obj, coeffs, T, rho)
        % Dynamic viscosity of ethanol, Sotiriadou et al. 2023   mu in (Pa-s)
        % -----------------------------------------------------------------
        % Arguments
        %   coeffs = 21 element row vector:
        %       (1:6)   b_i, i = 0 ... 5, numerator of eta0 in Tr
        %       (7:8)   c_i, i = 0 ... 1, denominator of eta0 in Tr
        %       (9)     epsilon/k_B (K)
        %       (10)    sigma (m)
        %       (11:17) d_i, i = 0 ... 6
        %       (18:19) d7, d8
        %       (20:21) e1, e2, coefficients of the residual term
        %   T = temperature (K)
        %   rho = density (kg/m^3)
        % -----------------------------------------------------------------
        % Assumptions
        %   1) The paper's viscosity is in uPa-s, converted to Pa-s here.
        %   2) The initial-density term eta1 is multiplied by rho in
        %       kg/m^3, as it was in regenEngineTradesCalculator.m
        % -----------------------------------------------------------------
        % Sources
        %   1) Sotiriadou, S., Ntonti, E., Velliadou, D., Antoniadis, K. D.,
        %       Assael, M. J., & Huber, M. L. (2023). Reference Correlation
        %       for the Viscosity of Ethanol from the Triple Point to 620 K
        %       and Pressures Up to 102 MPa
        % -----------------------------------------------------------------
            Tr = T ./ obj.T_crit;
            rho_r = rho ./ obj.rho_crit;

            % Zero-density viscosity
            numerator = zeros(size(Tr));
            for i = 1:6
                numerator = numerator + coeffs(i) .* Tr.^(i - 1);
            end
            denominator = zeros(size(Tr));
            for i = 1:2
                denominator = denominator + coeffs(6 + i) .* Tr.^(i - 1);
            end
            eta0 = numerator ./ denominator;

            % Second viscosity virial coefficient, reduced then dimensional
            Tstar = T ./ coeffs(9);
            Bstar = coeffs(18) .* Tstar.^-2.5 + coeffs(19) .* Tstar.^-5.5;
            for i = 0:6
                Bstar = Bstar + coeffs(11 + i) .* Tstar.^(-.25 .* i);
            end
            B_eta = Bstar .* Propellant.NA .* coeffs(10).^3 ./ obj.mlr_wgt;
            eta1 = eta0 .* B_eta;

            % Residual viscosity
            delta_eta = rho_r.^(2 ./ 3) .* Tr.^.5 .* (coeffs(20) .* rho_r + coeffs(21) .* (rho_r.^8 ./ (Tr.^4 .* (1 + rho_r.^2) - Tr.^2)));

            mu = (eta0 + eta1 .* rho + delta_eta) .* 1E-6;
        end

    end

    methods (Static, Access = private)

        function p_vap = modelAntoine(coeffs, T)
        % Vapor pressure by piecewise Antoine equation       p_vap in (Pa)
        % -----------------------------------------------------------------
        % Arguments
        %   coeffs = matrix, one row per valid temperature range:
        %       [T_low, T_high, A, B, C] (K, K, ~, K, K), range is
        %       T_low <= T < T_high. log10(p in bar) = A - B / (T + C)
        %   T = temperature (K)
        % -----------------------------------------------------------------
        % Assumptions
        %   1) T outside all ranges returns NaN, with a warning.
        % -----------------------------------------------------------------
        % Sources
        %   1) NIST Webbook
        % -----------------------------------------------------------------
            p_vap = nan(size(T));
            for i = 1:size(coeffs, 1)
                inRange = T >= coeffs(i, 1) & T < coeffs(i, 2);
                p_vap(inRange) = 10.^(coeffs(i, 3) - coeffs(i, 4) ./ (T(inRange) + coeffs(i, 5)));
            end
            % bar to Pa
            p_vap = p_vap .* 1E5;

            if any(T(:) < min(coeffs(:, 1)))
                warning("Some input T's are too low for the Antoine coefficients (below %g K)", min(coeffs(:, 1)));
            end
            if any(T(:) >= max(coeffs(:, 2)))
                warning("Some input T's are at or above the Antoine coefficients' upper limit (%g K)", max(coeffs(:, 2)));
            end
        end

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
                        error("Property ""%s"" is ""correlation"" in %s but ""%s"" and/or ""%s"" are missing", fieldName, Propellant.configFile, modelName, coeffsName);
                    end
                    data.(fieldName) = struct("model", data.(modelName), "coeffs", data.(coeffsName));
                    data = rmfield(data, [modelName, coeffsName]);
                end
            end
        end

        function listConfigNames()
        % Prints config names found in the config file
            configFileSearchMatches = dir(fullfile(pwd, "**", Propellant.configFile));
            if isempty(configFileSearchMatches)
                error("Config file ""%s"" not found in working directory ""%s""\n", Propellant.configFile, pwd);
            end
            rawText = readlines(fullfile(configFileSearchMatches(1).folder, configFileSearchMatches(1).name));
            allConfigsText = rawText(find(rawText == "%", 1) + 1:end);
            isHeaderLine = startsWith(allConfigsText, "\") & endsWith(allConfigsText, "\") & allConfigsText ~= "\END\";
            names = erase(allConfigsText(isHeaderLine), "\");
            fprintf("Valid propellant names in %s:\n", Propellant.configFile);
            fprintf("  %s\n", names);
        end

    end

end
