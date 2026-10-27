% The RocketEngine class object stores data about a physical rocket engine,
% like its geometry, material information, propellant information,
% performance metrics, etc. It takes design inputs, and computes combustion
% gas properties with RocketCEA (a Python wrapper for NASA CEA), as well as
% performance metrics, the physical engine contour (chamber, throat, and
% nozzle), coolant channel flow properties, wall temperatures, and stress
% estimates (physical and thermal).
% -------------------------------------------------------------------------
% Dependencies
%   1) Propellant.m
%   2) Material.m
%   3) RocketCEA 1.2.1 (Python library)
% -------------------------------------------------------------------------
% Assumptions
%   1) The engine uses 2 fluid propellants, one (1) oxidizer and one (1)
%       fuel
%   2) The engine uses regenerative cooling only. The fuel propellant is
%       used as the coolant. All of the coolant flows into a manifold inlet
%       near the engine exit and up channels in the engine sidewall to an
%       injector manifold interface, after which it enters the injector.
%   3) The coolant channel shape is square, and therefore the hydraulic
%       diameter of each channel is a side length of its square profile
%   4) Coolant properties are evaluated once, at the coolant inlet
%       temperature, and held constant along the channel length
%   5) Combustion gas flow is quasi 1D and isentropic. Gas specific heat
%       ratio and transport properties come from RocketCEA and are
%       evaluated at stations for a) the chamber, b) the throat, and c) the
%       nozzle
%   6) The physical engine contour is built from the following sections,
%       starting from the nozzle exit: a parabolic nozzle section (Rao
%       approximation), a circular diverging arc up to the throat, a
%       circular converging arc upstream of the throat, a parabolic
%       converging section in the chamber, then a straight chamber section
%       that makes up the rest of the chamber up to the injector faceplate
%       plane. The length of the chamber is calculated using L star. The
%       x-y plane is at the exit plane, and z increases in the direction
%       from the nozzle exit plane towards the injector faceplate plane
%   7) Heat transfer at each axial station is 1D (radial) and steady across
%       these media: 1) gas film (combusting flow), 2) hot wall
%       (cylindrical shell, engine side wall), and 3) coolant film (coolant
%       channel flow). The outer engine wall and the ambient are not
%       modeled; that is beyond the coolant is assumed adiabatic. Coolant
%       temperature rise along the channel is not modeled in the
%       station-wise analysis.
%   8) Hoop stress estimates use a thin wall assumption
%   9) Vectorization: all engines in an array are solved at the same time
%       and must share fuel, oxidizer and material. CEA is not
%       vectorizable,
%       so it is called once for each unique (pc, MR, expansion_ratio, CR).
% -------------------------------------------------------------------------
% Comments
%   1) All units are in SI, inputs and outputs
%   2) Constructing with array inputs makes an N-D array of engines. Each
%       input that has more than one element is a sweep dimension, in the
%       order of RocketEngine.sweepNames. Inputs with one element are held
%       constant. E.g. pc with 3 values and t_wall with 2 gives a 3x2
%       array.
%   3) Results are grouped in structs: gas (CEA), perf (performance), geom
%       (geometry and contour), coolant (inlet coolant state) and wall
%       (temperatures, heat flux, stresses). Station arrays are row
%       vectors, one element per axial station, from nozzle exit to
%       injector. Use extract() to pull one result out of an array of
%       engines as a numeric array.
%   4) Stations hold 8 B each and ~15 arrays per engine, so a stored array
%       of engines takes roughly 50 kB per engine.
% -------------------------------------------------------------------------
% Pending Additions/Changes
%   1) Amending assumption 7: Currently the temperate of the coolant at
%       each station is assumed to be the same, its value at the channel
%       inlet, and then the total coolant temperature rise is calculated
%       afterwards. This is not accurate, because at each station (going in
%       the +z direction) the coolant temperature would be higher, which
%       means the heat flux would be lower (same gas temperature, lower
%       temperature gradient = lower heat flux), which means less heat
%       would be getting pulled away from the hot wall. This means the code
%       underpredicts the hot wall temperature at stations upstream of the
%       throat with the error slightly growing as the stations approach the
%       injector faceplate plane. The change would be to change the code to
%       incrementally calculate each station's coolant flow temperature,
%       recalculate the heat flux, and subsequently the hot wall
%       temperature there.
%   2) Amending assumption 8: Code should be improved to use the thick wall
%       assumption because the geometry at the channels, with ribs on both
%       sides and a complicated inner/outer jacket structure is not as
%       simple as a plain thin walled pressure vessel. In fact, even a
%       thick wall assumption may not be sufficient, and a specific
%       derivation may be required.
%   3) Amending comment 1: Adding support for imperial units
%   4) Amending Nomenclature section: Once main first build of class
%       implementation is complete, fill out Nomenclature section
%   5) Amending GetAmbientPressure: Check that correlation is fitted
%       against tabulated values and do input checking to make sure
%       out-of-bounds values are rejected and an error is thrown
%   6) Amending comment 4: Recalculate this, maybe with some benchmark
%       examples
%   7) Amending constructor: Add RocketEngine_data.txt config file, and
%       make a default config for the constructor to pull data from
%   8) Amending constructor: Make impossible-input handling cases, so if
%       the inputs vary over ranges and some of those combinations produce
%       impossible engines, just create a log file for those and skip the
%       rest of their computations. Additionally, expand which variables
%       the user can sweep through
% -------------------------------------------------------------------------
% Nomenclature
%   pc = chamber pressure (Pa)
%   MR = oxidizer to fuel mixture ratio by mass
%   Tr = reduced temperature
%   T_wg = hot wall temperature at the gas side (K)
%   T_wl = hot wall temperature at the coolant side (K)
%   q = heat flux (W/m^2)
%   <Symbol> = <Definition> <(Units)>
% -------------------------------------------------------------------------
% Document Version 1.0
% -------------------------------------------------------------------------
% MATLAB Version R2024b
% -------------------------------------------------------------------------
% Developed by Alex Vance (AlexVance00 on Github)
classdef RocketEngine

    %% Properties
    properties (Constant, Hidden)
        % Physical constants
        g0 = 9.81; % Gravitational acceleration on Earth at sea level, m/s^2
        R_universal = 8.314; % Universal gas constant, J/kg-K
        p_SL = 101325; % Total air pressure of still air at sea level, Pa
        T_SL = 288.15; % Total temperature of still air at sea level, K
        mlr_wgt_air = 0.02896968; % Molar weight of air, kg/mol
        lapse_rate = 0.00976; % Temperature lapse rate with altitude used for ambient pressure calculation, K/m

        % Number of axial stations in each section of the contour. These set
        % the size of every station array, so they're shared by all engines
        numStations = struct("nozzleParabolic", 5, ...
            "nozzleCircular", 5, "chamberCircular", 5, ...
            "chamberParabolic", 5, "chamberLinear", 5);

        % Reynolds number at which pipe flow is considered to transition
        % from laminar to turbulent with the pipe diameter as the
        % characteristic length
        Re_D_transition = 2300;
        % Stations whose wall slope |dr/dz| is below this are treated as
        % straight, so the radius of curvature term in Bartz drops out
        tol_R_curv = 0.05;
        % Stations whose area ratio is within this of 1 are at the throat
        tol_throat = 0.0001;
        % Wall temperature convergence tolerance (K) and max iterations
        tol_T_wg = 5;
        maxIter_T_wg = 20;

        % Names of the inputs that can be swept, in the order of the dimensions
        % of the engine array
        sweepNames = ["pc", "MR", "expansion_ratio", "mdot", "d_channel", ...
            "num_channels", "T_coolant_i", "t_wall", "CR"];
    end

    properties (SetAccess = private)
        % Identity
        name = "Hydra";

        % Environment
        alt = 0; % Altitude above sea level
        T_amb = RocketEngine.T_SL; % Ambient temperature, K
        p_amb = RocketEngine.p_SL; % Ambient pressure, Pa

        % Information about the material the engine is made of
        material Material % Material class object

        % Information about the engine propellants
        fuel Propellant % Engine fuel propellant, Propellant class object
        oxidizer Propellant % Engine oxidizer propellant, Propellant class object

        % Design inputs
        pc = []; % Chamber pressure, Pa
        MR = []; % Oxidizer to fuel mass ratio
        expansion_ratio = []; % Nozzle exit area to throat area
        mdot = []; % Total propellant mass flow rate, kg/s
        d_channel = []; % side length of coolant channel square profile, m
        num_channels = []; % Number of coolant channels
        T_coolant_i = []; % Coolant temperature at channel inlet, K
        t_wall = []; % Hot wall thickness, m

        % Parameters for combustion chamber design
        CR = []; % Contraction ratio, chamber area / throat area
        Lstar = []; % Characteristic chamber length (proxy for residence time, propellant-dependent), chamber volume / throat area, m
        etaRc = .2; % How far (fraction) between throat and chamber radius the converging arc ends

        % Design factors
        DF_stress = 2;
        DF_temp = 1.5; % Design factor for gas-side hot wall temperature
        DF_p_coolant = 1.2; % Design factor for minimum coolant pressure to prevent boiling in coolant channels

        % Parameters for Rao approximation for bell nozzle contour design
        fpct = 1; % Nozzle length as fraction of a 15 deg half-angle conical nozzle
        thetaN = 20; % Nozzle parabola initial angle used in , deg

        % Results, see Comments in the class header
        gasProps struct % Combustion gas state from CEA
        perf struct % Performance
        geom struct % Geometry and contour
        coolant struct % Coolant state at the channel inlet
        wall struct % Wall temperatures, heat flux and stresses

        % Flags
        flagIsSolved logical = false; % Whether solve() has run for this engine
    end

    %% Public methods
    methods (Access = public)

        % Constructor
        function obj = RocketEngine(varargin)
        % Makes a RocketEngine, or an N-D array of them if any design input
        % has more than one element, and solves it unless told not to
        %                                        obj, array of RocketEngine
        % -----------------------------------------------------------------
        % Arguments (all optional, called as RocketEngine(pc = ..., ...))
        %   fuel, oxidizer = Propellant config names (string) or Propellant
        %       objects
        %   material = Material config name (string) or Material object
        %   alt = altitude (m), T_amb = ambient temperature (K),
        %   p_amb = ambient pressure (Pa), empty to compute from alt
        %   pc, MR, expansion_ratio, mdot, d_channel, num_channels,
        %   T_coolant_i, t_wall, CR, Lstar = see properties.
        %       T_coolant_i empty -> T_amb.
        %       Any of these may be an array to make an array of engines.
        %   thetaN, fpct, etaRc, DF_stress, DF_temp, DF_p_coolant = see
        %       properties, scalars
        %   solve = false to only make the objects without running the
        %       analysis (default true)            
            args = RocketEngine.parseArgs(varargin{:});

            % No args returns empty object
            if nargin == 0
                return;
            end

            % Propellants and material, from config name or passed in
            % object
            fuel = RocketEngine.asObject(args.fuel, "Propellant");
            oxidizer = RocketEngine.asObject(args.oxidizer, "Propellant");
            material = RocketEngine.asObject(args.material, "Material");
            if fuel.type ~= "FUEL"
                error("Propellant ""%s"" is type ""%s"", but was passed as the fuel", fuel.name, fuel.type);
            end
            if oxidizer.type ~= "OXIDIZER"
                error("Propellant ""%s"" is type ""%s"", but was passed as the oxidizer", oxidizer.name, oxidizer.type);
            end

            % Inputs that default to something that depends on other inputs
            if isnan(args.p_amb)
                args.p_amb = RocketEngine.GetAmbientPressure(alt = args.alt);
            end

            % Sweep: every input in sweepNames with more than one element
            % is a dimension of the array of engines
            sweepNames = RocketEngine.sweepNames;
            values = cell(1, numel(sweepNames));
            for j = 1:numel(sweepNames)
                values{j} = args.(sweepNames(j))(:);
                if isempty(values{j})
                    error("Input ""%s"" can't be empty", sweepNames(j));
                end
            end
            sweepIndices = find(cellfun(@numel, values) > 1);
            if isempty(sweepIndices)
                gridSize = [1, 1];
                grids = {};
            else
                grids = cell(1, numel(sweepIndices));
                [grids{:}] = ndgrid(values{sweepIndices});
                gridSize = size(grids{1});
            end
            numEngines = prod(gridSize);

            % Preallocate array of obj types for speed's sake, 1 element
            % for each unique combination of input parameters
            if numEngines ~= 1
                obj(numEngines) = RocketEngine();
            end

            % Assign args
            for i = 1:numEngines
                obj(i).name = args.name;
                obj(i).alt = args.alt;
                obj(i).T_amb = args.T_amb;
                obj(i).p_amb = args.p_amb;
                obj(i).material = material;
                obj(i).fuel = fuel;
                obj(i).oxidizer = oxidizer;
                for j = 1:numel(sweepNames)
                    gridIndex = find(sweepIndices == j);
                    if isempty(gridIndex)
                        obj(i).(sweepNames(j)) = values{j};
                    else
                        obj(i).(sweepNames(j)) = grids{gridIndex}(i);
                    end
                end
                obj(i).Lstar = args.Lstar;
                obj(i).etaRc = args.etaRc;
                obj(i).DF_stress = args.DF_stress;
                obj(i).DF_temp = args.DF_temp;
                obj(i).DF_p_coolant = args.DF_p_coolant;
                obj(i).fpct = args.fpct;
                obj(i).thetaN = args.thetaN;
            end

            if numEngines ~= 1
                obj = reshape(obj, gridSize);
            end

            if ~args.flagIsSolved
                obj = obj.solve();
            end

            return;
        end

        function objs = solve(objs)
        % Runs the whole analysis for every engine in the array, all at
        % once
        %                                   objs, array of RocketEngine
            inputs = RocketEngine.organizeInputs(objs);

            % From assumption that only one unique fuel, oxidizer, and
            % engine material are being run
            fuel = objs(1).fuel;
            oxidizer = objs(1).oxidizer;
            material = objs(1).material;
            numEngines = numel(objs);

            %% Combustion gas
            gasProps = RocketEngine.runCEA(inputs, fuel, oxidizer);

            %% Performance and sizing
            ge = gasProps.gamma_e;
            perf.pe = gasProps.pe;
            perf.Cf = sqrt(2 .* ge.^2 ./ (ge - 1) .* (2 ./ (ge + 1)).^((ge + 1) ./ (ge - 1)) .* (1 - (gasProps.pe ./ inputs.pc).^((ge - 1) ./ ge))) + inputs.expansion_ratio .* ((gasProps.pe - inputs.p_amb) ./ inputs.pc);
            perf.flow_sonic = sqrt(gasProps.gamma_t .* RocketEngine.R_universal ./ gasProps.MW_t .* gasProps.T_t);
            perf.At = inputs.mdot ./ perf.flow_sonic ./ gasProps.rho_t;
            perf.dt = sqrt(perf.At .* 4 ./ pi);
            perf.Ae = perf.At .* inputs.expansion_ratio;
            perf.de = sqrt(perf.Ae .* 4 ./ pi);
            perf.Vc = perf.At .* inputs.Lstar;
            perf.thrust = perf.Cf .* perf.At .* inputs.pc;
            perf.prop_cost_rate = inputs.mdot .* (oxidizer.cost .* inputs.MR ./ (1 + inputs.MR) + fuel.cost ./ (1 + inputs.MR));
            perf.Isp = perf.thrust ./ inputs.mdot ./ RocketEngine.g0;
            perf.cstar = perf.Isp .* RocketEngine.g0 ./ perf.Cf;
            perf.eta_cstar = perf.cstar ./ gasProps.cstar_theo;

            %% Contour
            geom = RocketEngine.buildContour(inputs, perf.dt, perf.Vc);
            geom.vol_engine = RocketEngine.getShellVolume(geom, inputs);
            perf.TWR = perf.thrust ./ (RocketEngine.g0 .* material.rho .* geom.vol_engine);

            %% Combustion gas state along the contour
            gasLine = RocketEngine.getGasAlongContour(geom, gasProps);

            %% Coolant at the channel inlet
            coolant = RocketEngine.getCoolantInlet(inputs, fuel);

            %% Throat estimate (lumped, Dittus-Boelter on both sides)
            wall = RocketEngine.getThroatEstimate(inputs, gasProps, perf, coolant);
            % Crude worst case coolant temperature rise: throat heat flux
            % at every cross section of the whole channel
            % (kept as it was in regenEngineTradesCalculator.m)
            l_channel = 2 .* max(geom.z, [], 2);
            deltaT_coolantslice = wall.q_throat .* inputs.d_channel ./ coolant.mdot_channel ./ coolant.cp;
            wall.T_coolant_f = inputs.T_coolant_i + deltaT_coolantslice .* l_channel;
            % Minimum coolant pressure to keep it from boiling in the channels
            wall.P_coolant_min = fuel.getPvap(wall.T_coolant_f) ./ inputs.DF_p_coolant;

            %% Station-wise wall temperatures and heat flux
            stations = RocketEngine.solveWallTemperatures(inputs, material, gasProps, perf, geom, gasLine, coolant);
            wall = RocketEngine.mergeStructs(wall, stations);

            %% Stress and design limits
            wall.hoop_stress_doghouse = inputs.pc .* geom.Rc ./ inputs.t_wall;
            wall.hoop_stress_fins = inputs.pc .* (geom.Rc + inputs.t_wall) ./ inputs.d_channel;
            wall.therm_stress = material.E .* material.CTE .* (wall.T_wg - wall.T_wl);
            wall.total_stress_doghouse = wall.therm_stress + wall.hoop_stress_doghouse;
            wall.total_stress_fins = wall.therm_stress + wall.hoop_stress_fins;
            wall.yield_stress = material.getSigmaY(wall.T_wg);
            wall.yield_stress_allowable = wall.yield_stress ./ inputs.DF_stress;
            % Design factor applied to the melting temperature in deg C, as
            % regenEngineTradesCalculator.m did, then back to K
            wall.T_wall_allowable = (material.T_melt - 273.15) ./ inputs.DF_temp + 273.15;

            %% Store results in the engines
            gasLineFields = struct("M", gasLine.M, "T_gas", gasLine.T);
            geom = RocketEngine.mergeStructs(geom, gasLineFields);
            for i = 1:numEngines
                objs(i).gasProps = RocketEngine.getRow(gasProps, i);
                objs(i).perf = RocketEngine.getRow(perf, i);
                objs(i).geom = RocketEngine.getRow(geom, i);
                objs(i).coolant = RocketEngine.getRow(coolant, i);
                objs(i).wall = RocketEngine.getRow(wall, i);
                objs(i).flagIsSolved = true;
            end

            return;
        end

        function values = extract(objs, resultPath)
        % Pulls one result out of every engine in an array
        %   values, array the same size as objs for a scalar result, or the
        %   size of objs with one more trailing dimension (stations) for a
        %   station array result
        % -----------------------------------------------------------------
        % Arguments
        %   resultPath = name of the property, or of the result struct
        %       followed by a dot and a field, e.g. "pc" or "perf.thrust" or
        %       "wall.T_wg"
        % -----------------------------------------------------------------
            arguments
                objs RocketEngine
                resultPath (1, 1) string
            end

            pathParts = split(resultPath, ".");
            cellValues = cell(size(objs));
            for i = 1:numel(objs)
                value = objs(i);
                for j = 1:numel(pathParts)
                    value = value.(pathParts(j));
                end
                cellValues{i} = value;
            end

            numPerEngine = numel(cellValues{1});
            stacked = reshape(vertcat(cellValues{:}), numel(objs), numPerEngine);
            if numPerEngine == 1
                values = reshape(stacked, size(objs));
            else
                values = reshape(stacked, [size(objs), numPerEngine]);
            end
        end

    end

    %% Physics: fluid dynamics and heat transfer (folded in from the Get*.m files)
    methods (Static, Access = public)

        function printValidMaterialNames()
        % Prints the names of every material that can be passed as the
        % material input, the configs. in material_data.txt
            Material.printValidNames();
        end

        function printValidPropellantNames()
        % Prints the names of every propellant that can be passed as the
        % fuel or oxidizer input, the configs. in propellant_data.txt
            Propellant.printValidNames();
        end

        function p_amb = GetAmbientPressure(args)
        % Gets ambient pressure at altitude
        % -----------------------------------------------------------------
        % Arguments
        %   alt = altitude above sea level (m)
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Temperature falls linearly with altitude at rate of
        %       lapse_rate
        % -----------------------------------------------------------------
            arguments
                args.alt = [];
            end
            RocketEngine.requireArgs(args, "alt");
            alt = args.alt;

            p_amb = RocketEngine.p_SL .* (1 + (RocketEngine.lapse_rate .* alt) ./ RocketEngine.T_SL).^(-RocketEngine.g0 .* RocketEngine.mlr_wgt_air ./ (RocketEngine.R_universal .* RocketEngine.lapse_rate));
        end

        function sigma = GetBartzCorrectionFactor(args)
        % Correction factor "sigma" for the GetBartzh_g method, for property
        % variations across the boundary layer                  sigma in (~)
        % -----------------------------------------------------------------
        % Arguments
        %   T_wg = gas-side hot wall temperature (units must match T_c)
        %   T_c = chamber stagnation temperature (units must match T_wg)
        %   gammay = specific heat ratio
        %   M = flow Mach number
        %   All the same size, or scalars
        % -----------------------------------------------------------------
        % Assumptions
        %   See Assumption 4 of the class header
        % -----------------------------------------------------------------
        % Sources
        %   1) Design of Liquid Propellant Rocket Engines 2nd. Ed. by Huzel
        %       and Huang p. 101 Eq. (4-14) (Bartz Correlation)
        % -----------------------------------------------------------------
            arguments
                args.T_wg = [];
                args.T_c = [];
                args.gammay = [];
                args.M = [];
            end
            RocketEngine.requireArgs(args, ["T_wg", "T_c", "gammay", "M"]);
            T_wg = args.T_wg;
            T_c = args.T_c;
            gammay = args.gammay;
            M = args.M;

            sigma = 1 ./ ...
                (((1 ./ 2) .* (T_wg ./ T_c) .* ...
                    (1 + (gammay - 1) ./ 2 .* M.^2) + 1 ./ 2).^.68 .* ...
                (1 + (gammay - 1) ./ 2 .* M.^2).^.12);
        end

        function h_g = GetBartzh_g(args)
        % Gas-side convective heat transfer coefficient at the hot wall
        % by the Bartz correlation                        h_g in (W/m^2-K)
        % -----------------------------------------------------------------
        % Arguments
        %   D_t = throat diameter (m)
        %   mu = gas viscosity (kg/m-s)
        %   c_p = gas specific heat at constant pressure (J/kg-K)
        %   Pr = gas Prandtl number
        %   p_c = chamber pressure (Pa)
        %   g = gravitational acceleration (m/s^2)
        %   c_star = characteristic velocity (m/s)
        %   R = radius of curvature of the wall at this station (m)
        %   A_t = throat area (units must match A)
        %   A = area at this station (units must match A_t)
        %   sigma = correction factor, from GetBartzCorrectionFactor
        %   All the same size, or scalars
        % -----------------------------------------------------------------
        % Assumptions
        %   1) The correlation is in English units, so inputs are converted
        %       from SI to English inside, and the result back to SI
        %   See also Assumption 4 of the class header
        % -----------------------------------------------------------------
        % Sources
        %   1) Design of Liquid Propellant Rocket Engines 2nd. Ed. by Huzel
        %       and Huang p. 100 Eq. (4-13) (Bartz Correlation)
        % -----------------------------------------------------------------
            arguments
                args.D_t = [];
                args.mu = [];
                args.c_p = [];
                args.Pr = [];
                args.p_c = [];
                args.g = [];
                args.c_star = [];
                args.R = [];
                args.A_t = [];
                args.A = [];
                args.sigma = [];
            end
            RocketEngine.requireArgs(args, ["D_t", "mu", "c_p", "Pr", "p_c", "g", "c_star", "R", "A_t", "A", "sigma"]);
            Pr = args.Pr;
            A_t = args.A_t;
            A = args.A;
            sigma = args.sigma;

            % Unit conversions
            % m to in
            D_t = args.D_t .* 39.3701;
            R = args.R .* 39.3701;
            % kg/m-sec to lbm/in-sec
            mu = args.mu .* 0.0254 ./ 0.45359237;
            % Pa to psi
            p_c = args.p_c ./ 6894.76;
            % m/sec to ft/sec
            c_star = args.c_star .* 3.28084;
            % m/sec^2 to ft/sec^2
            g = args.g .* 3.28084;
            % J/kg-K to Btu/lbm-deg F
            c_p = args.c_p .* 0.000947817 ./ 2.20462 ./ 1.8;

            % Bartz correlation (yields h_g in Btu/in^2-sec-deg F)
            h_g = ((.026 ./ D_t.^.2) .* (mu.^.2 .* c_p ./ Pr.^.6) .* (p_c .* g ./ c_star).^.8 .* (D_t ./ R).^.1) .* (A_t ./ A).^.9 .* sigma;

            % Btu/in^2-sec-deg F to J/m^2-sec-K
            h_g = h_g .* 1055.06 ./ 0.0254.^2 ./ (5 ./ 9);
        end

        function f = GetDarcyFrictionFactor(args)
        % Darcy friction factor from a pressure drop along a channel   f in (~)
        % -----------------------------------------------------------------
        % Arguments
        %   D = hydraulic diameter of channel (units must match L)
        %   L = channel length (units must match D)
        %   delta_p = pressure change from inlet to outlet (Pa)
        %   rho = channel flow density (kg/m^3)
        %   u_bar = channel flow average speed (m/s)
        % -----------------------------------------------------------------
        % Assumptions
        %   See Assumptions 2 and 3 of the class header
        % -----------------------------------------------------------------
            arguments
                args.D = [];
                args.L = [];
                args.delta_p = [];
                args.rho = [];
                args.u_bar = [];
            end
            RocketEngine.requireArgs(args, ["D", "L", "delta_p", "rho", "u_bar"]);

            f = args.D ./ args.L .* args.delta_p ./ (args.rho .* args.u_bar.^2 ./ 2);
        end

        function f = GetHaalandFrictionFactor(args)
        % Darcy friction factor by the Haaland equation                f in (~)
        % -----------------------------------------------------------------
        % Arguments
        %   epsilon = surface roughness (m)
        %   D = hydraulic diameter of channel (m)
        %   Re = channel flow Reynolds number
        % -----------------------------------------------------------------
        % Assumptions
        %   See Assumptions 2 and 3 of the class header
        % -----------------------------------------------------------------
            arguments
                args.epsilon = [];
                args.D = [];
                args.Re = [];
            end
            RocketEngine.requireArgs(args, ["epsilon", "D", "Re"]);

            f = 1 ./ (-1.8 .* log10((args.epsilon ./ args.D ./ 3.7).^1.11 + 6.9 ./ args.Re)).^2;
        end

        function x_fd_h_over_D = GetHydrodynamicEntryLength(args)
        % Internal flow hydrodynamic entry length, nondimensionalized by the
        % tube diameter                                     x_fd_h_over_D in (~)
        % -----------------------------------------------------------------
        % Arguments
        %   Re_D = Reynolds number using tube diameter, any size
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Laminar below Re_D_transition, turbulent above
        %   in x_fd_h_over_D: x = axial distance from datum, fd = fully
        %       developed, h = hydrodynamic, D = diameter of tube
        % -----------------------------------------------------------------
        % Sources
        %   1) Laminar: Introduction to Heat Transfer 3rd Ed. by Incropera
        %       and DeWitt, p. 389, Eq. (8.3)
        %   2) Turbulent: Fundamentals of Thermal Fluid Sciences by Cengel
        %       and Turner, p. 760, Eq. (17-42), rearranged to leave out D
        % -----------------------------------------------------------------
            arguments
                args.Re_D = [];
            end
            RocketEngine.requireArgs(args, "Re_D");
            Re_D = args.Re_D;

            x_fd_h_over_D = 4.4 .* Re_D.^(1 ./ 6);
            laminar = Re_D <= RocketEngine.Re_D_transition;
            x_fd_h_over_D(laminar) = .05 .* Re_D(laminar);
        end

        function x_fd_t_over_D = GetThermalEntryLength(args)
        % Internal flow thermal entry length, nondimensionalized by the tube
        % diameter                                          x_fd_t_over_D in (~)
        % -----------------------------------------------------------------
        % Arguments
        %   Re_D = Reynolds number using tube diameter, any size
        %   Pr = Prandtl number, same size as Re_D, or scalar
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Laminar below Re_D_transition, turbulent above
        %   2) The turbulent expression doesn't depend on Pr or D. The
        %       earlier GetThermalEntryLength.m multiplied it by D, which
        %       is a units mistake for a nondimensional length (the comment
        %       there said it was meant to be left out), so it's left out here
        %   in x_fd_t_over_D: x = axial distance from datum, fd = fully
        %       developed, t = thermal, D = diameter of tube
        % -----------------------------------------------------------------
        % Sources
        %   1) Laminar: Introduction to Heat Transfer 3rd Ed. by Incropera
        %       and DeWitt, p. 393, Eq. (8.23)
        %   2) Turbulent: Fundamentals of Thermal Fluid Sciences by Cengel
        %       and Turner, p. 760, Eq. (17-42), rearranged to leave out D
        % -----------------------------------------------------------------
            arguments
                args.Re_D = [];
                args.Pr = [];
            end
            RocketEngine.requireArgs(args, ["Re_D", "Pr"]);
            Re_D = args.Re_D;
            Pr = args.Pr .* ones(size(Re_D));

            x_fd_t_over_D = 4.4 .* Re_D.^(1 ./ 6);
            laminar = Re_D <= RocketEngine.Re_D_transition;
            x_fd_t_over_D(laminar) = .05 .* Re_D(laminar) .* Pr(laminar);
        end

        function Nu_D = GetGnielinskiNu_D(args)
        % Nusselt number (tube diameter as characteristic length) by the
        % Gnielinski correlation                                  Nu_D in (~)
        % -----------------------------------------------------------------
        % Arguments
        %   f = friction factor
        %   Re_D = Reynolds number using channel diameter
        %   Pr = Prandtl number
        %   D = channel hydraulic diameter (m)
        %   L = channel length (m)
        %   x = position in channel relative to its start (m), any size,
        %       the others are the size of x or compatible with it
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Turbulent flow, Re_D < 5E6, 0 < D/L < 1, properties at T_m =
        %       (T_1 + T_2)/2
        %   2) Stations that are still in the entrance region (x/D below the
        %       longer of the thermal and hydrodynamic entry lengths) get the
        %       entrance correction (1 + (D/L)^(2/3)), the others are
        %       fully developed
        %   3) L = 0 would divide by zero, so it's replaced by 2.5% of D
        % -----------------------------------------------------------------
        % Sources
        %   1) Heat Convection 2nd Ed. by Latif M. Jiji, p. 404, Eq. (10.17)
        %       (Gnielinski correlation)
        % -----------------------------------------------------------------
            arguments
                args.f = [];
                args.Re_D = [];
                args.Pr = [];
                args.D = [];
                args.L = [];
                args.x = [];
            end
            RocketEngine.requireArgs(args, ["f", "Re_D", "Pr", "D", "L", "x"]);
            f = args.f;
            Re_D = args.Re_D;
            Pr = args.Pr;
            D = args.D;
            x = args.x;

            % Avoid divide by zero if L = 0. Scale to D
            L = args.L + (args.L == 0) .* (.01 .* 2.5 .* D);

            % Check if developing, using entrance lengths
            x_fd_t_over_D = RocketEngine.GetThermalEntryLength(Re_D = Re_D, Pr = Pr);
            x_fd_h_over_D = RocketEngine.GetHydrodynamicEntryLength(Re_D = Re_D);
            developing = (x ./ D) < max(x_fd_t_over_D, x_fd_h_over_D);

            D_over_L = developing .* (D ./ L);

            Nu_D = ((f ./ 8) .* (Re_D - 1000) .* Pr) ./ (1 + 12.7 .* ...
                sqrt((f ./ 8)) .* (Pr.^(2 ./ 3) - 1)) .* (1 + D_over_L.^(2 ./ 3));
        end

        function h = Geth(args)
        % Convective heat transfer coefficient from a Nusselt number
        %                                                  h in (W/m^2-K)
        % -----------------------------------------------------------------
        % Arguments
        %   Nu_D = Nusselt number using tube diameter as the characteristic
        %       length
        %   D = channel hydraulic diameter (m)
        %   k = fluid thermal conductivity (W/m-K)
        % -----------------------------------------------------------------
            arguments
                args.Nu_D = [];
                args.k = [];
                args.D = [];
            end
            RocketEngine.requireArgs(args, ["Nu_D", "k", "D"]);

            h = args.Nu_D .* args.k ./ args.D;
        end

        function R = GetRadiusOfCurvature(args)
        % Radius of curvature of the wall contour from finite differences
        %                                                           R in (m)
        % -----------------------------------------------------------------
        % Arguments
        %   dx = finite delta x (units must match dy)
        %   dy = finite delta y (units must match dx)
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Make sure dx, dy orders of magnitude are >= 2 lower than the
        %       chamber radius
        %   2) In the engine, dx is delta r and dy is delta z. Note that the
        %       formula takes (dy/dx)^2 as the second derivative, as
        %       GetRadiusOfCurvature.m did
        % -----------------------------------------------------------------
            arguments
                args.dx = [];
                args.dy = [];
            end
            RocketEngine.requireArgs(args, ["dx", "dy"]);
            dx = args.dx;
            dy = args.dy;

            R = (1 + (dy ./ dx).^2).^(3 ./ 2) ./ abs(dy.^2 ./ dx.^2);
        end

        function TR = GetThermalResistance(args)
        % Thermal resistance for one layer of the wall stack, per unit area
        % at a reference radius                                TR in (m^2-K/W)
        % -----------------------------------------------------------------
        % Arguments
        %   heat_transfer_type = "conduction" or "forced convection"
        %       ("free convection" isn't implemented)
        %   k = thermal conductivity (W/m-K), for conduction
        %   h = convective heat transfer coefficient (W/m^2-K), for
        %       convection
        %   r_0 = reference radius for which the basis heat flux is defined
        %   r_1 = for conduction: starting radius, for convection: radius
        %       of the boundary
        %   r_2 = ending radius, for conduction
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Accounts for the unit-width increase of an angular cross
        %       section with radius, normalized to the reference radius
        %       (see Assumption 6 of the class header)
        % -----------------------------------------------------------------
            arguments
                args.heat_transfer_type (1, 1) string = ""
                args.k = [];
                args.h = [];
                args.r_0 = [];
                args.r_1 = [];
                args.r_2 = [];
            end

            switch args.heat_transfer_type
                case "conduction"
                    RocketEngine.requireArgs(args, ["k", "r_0", "r_1", "r_2"]);
                    TR = 1 ./ args.k .* args.r_0 .* log(args.r_2 ./ args.r_1);
                case "forced convection"
                    RocketEngine.requireArgs(args, ["h", "r_0", "r_1"]);
                    TR = 1 ./ args.h .* args.r_0 ./ args.r_1;
                case "free convection"
                    error("'%s' value for 'heat_transfer_type' parameter not yet implemented", args.heat_transfer_type);
                otherwise
                    error("Invalid input for 'heat_transfer_type' parameter");
            end
        end

        function M = GetMachFromAreaRatio(args)
        % Mach number from the isentropic area ratio, by Newton's method,
        % vectorized (replaces the symbolic Newton loop in the old calculator)
        %                                                             M in (~)
        % -----------------------------------------------------------------
        % Arguments
        %   area_ratio = A / A_throat, any size
        %   gamma = specific heat ratio, same size as area_ratio, or scalar
        %   subsonic = true where the subsonic root is wanted, false where the
        %       supersonic one is, same size as area_ratio, or scalar
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Isentropic, 1D, calorically perfect gas flow (Assumption 4 of
        %       the class header)
        %   2) Stations with area ratio within tol_throat of 1 return M = 1
        % -----------------------------------------------------------------
            arguments
                args.area_ratio = [];
                args.gamma = [];
                args.subsonic = [];
            end
            RocketEngine.requireArgs(args, ["area_ratio", "gamma", "subsonic"]);
            areaRatio = args.area_ratio;
            gammay = args.gamma .* ones(size(areaRatio));
            subsonic = args.subsonic & true(size(areaRatio));

            M = ones(size(areaRatio));
            free = abs(areaRatio - 1) > RocketEngine.tol_throat;
            ARv = areaRatio(free);
            gv = gammay(free);

            % Starting guesses
            Mv = 1.2 .* ones(size(ARv));
            Mv(subsonic(free)) = .01;

            notDone = false(size(ARv));
            for iter = 1:200
                a = (gv + 1) ./ (2 .* (gv - 1));
                B = 2 ./ (gv + 1) .* (1 + (gv - 1) ./ 2 .* Mv.^2);
                residual = B.^a ./ Mv - ARv;
                notDone = abs(residual) >= 1E-9;
                if ~any(notDone)
                    break;
                end
                dB = 2 ./ (gv + 1) .* (gv - 1) .* Mv;
                slope = a .* B.^(a - 1) .* dB ./ Mv - B.^a ./ Mv.^2;
                Mv(notDone) = max(Mv(notDone) - residual(notDone) ./ slope(notDone), 1E-6);
            end
            if any(notDone)
                warning("Mach number from area ratio didn't converge at %d stations", nnz(notDone));
            end
            M(free) = Mv;
        end

    end

    %% Stages of solve()
    methods (Static, Access = private)

        function args = parseArgs(args)
                arguments
                    args.name = "Hydra"
                    args.fuel = "ETHANOL_L"
                    args.oxidizer = "N2O_L"
                    args.material = "AlSi10Mg"
                    args.alt (1, 1) double = 342
                    args.T_amb (1, 1) double = 273.15 + 20
                    args.p_amb (1, 1) double = NaN
                    args.pc double = 128 * 6894.76
                    args.MR double = 1.6
                    args.expansion_ratio double = 1.8
                    args.mdot double = 1.2 * .45359237
                    args.d_channel double = 1.2E-3
                    args.num_channels double = 48
                    args.T_coolant_i double = []
                    args.t_wall double = .8E-3
                    args.CR double = 6
                    args.Lstar double = .6
                    args.thetaN (1, 1) double = 20
                    args.fpct (1, 1) double = 1
                    args.etaRc (1, 1) double = .2
                    args.DF_stress (1, 1) double = 2
                    args.DF_temp (1, 1) double = 1.5
                    args.DF_p_coolant (1, 1) double = 1.2
                    args.flagIsSolved (1, 1) logical = false
                end
            end

        function requireArgs(args, requiredNames)
        % Throws error if a required name-value argument was left unassigned
            for requiredName = requiredNames
                if isempty(args.(requiredName))
                    error("No input for non-optional '%s' argument", requiredName);
                end
            end
        end

        function obj = asObject(value, className)
        % Makes a Propellant or Material from a config name, or passes
        % through one that's already made
            if isa(value, className)
                obj = value;
            elseif isstring(value) || ischar(value)
                obj = feval(className, string(value));
            else
                error("Expected a %s config name or object, but got a %s object of size [%s] starting with %d", className, class(value), join(string(size(value)), ", "), value(1));
            end
            if ~isscalar(obj)
                error("Expected a single %s", className);
            end
        end

        function inputs = organizeInputs(objs)
        % Collects every engine's design inputs into column vectors
            objs = objs(:);
            fuelNames = arrayfun(@(e) char(e.fuel.name), objs, UniformOutput = false);
            oxidizerNames = arrayfun(@(e) char(e.oxidizer.name), objs, UniformOutput = false);
            materialNames = arrayfun(@(e) char(e.material.name), objs, UniformOutput = false);
            if numel(unique(fuelNames)) > 1 || numel(unique(oxidizerNames)) > 1 || numel(unique(materialNames)) > 1
                error("All engines solved together must have the same fuel, oxidizer and material");
            end

            inputNames = ["alt", "T_amb", "p_amb", "pc", "MR", "expansion_ratio", ...
                "mdot", "d_channel", "num_channels", "T_coolant_i", "t_wall", "CR", "Lstar", "thetaN", "fpct", "etaRc", ...
                "DF_stress", "DF_temp", "DF_p_coolant"];
            for inputName = inputNames
                inputs.(inputName) = reshape([objs.(inputName)], [], 1);
            end
        end

        function gasProps = runCEA(inputs, fuel, oxidizer)
        % Gets combustion gas state from CEA for every engine, calling CEA
        % once per unique (pc, MR, expansion_ratio, CR)
        % -----------------------------------------------------------------
        % Assumptions
        %   1) CEA's inputs and outputs are English units, converted here
        % -----------------------------------------------------------------
            % Import RocketCEA library
            py.importlib.import_module("rocketcea");

            CEAkeys = [inputs.pc, inputs.MR, inputs.expansion_ratio, inputs.CR];
            [uniqueCEAKeys, ~, mapToUnique] = unique(CEAkeys, "rows");
            numUniqueCEARuns = size(uniqueCEAKeys, 1);

            fieldNames = ["MW_t", "gamma_t", "gamma_e", "gamma_c", "pe", ...
                "Isp_theo", "T_c", "T_t", "T_e", "rho_t", "cp_t", "mu_t", ...
                "k_t", "Pr_t", "cp_c", "mu_c", "k_c", "Pr_c", "cstar_theo"];
            uniqueCEARuns = struct();
            for fieldName = fieldNames
                uniqueCEARuns.(fieldName) = zeros(numUniqueCEARuns, 1);
            end

            % 
            ceaByCR = containers.Map("KeyType", "double", "ValueType", "any");
            for u = 1:numUniqueCEARuns
                pc = uniqueCEAKeys(u, 1);
                MR = uniqueCEAKeys(u, 2);
                expansion_ratio = uniqueCEAKeys(u, 3);
                CR = uniqueCEAKeys(u, 4);
                if ~isKey(ceaByCR, CR)
                    ceaByCR(CR) = py.rocketcea.cea_obj.CEA_Obj(fuelName = fuel.name, oxName = oxidizer.name, fac_CR = CR);
                end
                cea = ceaByCR(CR);
                % Pa to psi
                pc = pc * 0.000145038;

                throatProps = double(cea.get_Throat_MolWt_gamma(Pc = pc, MR = MR, eps = expansion_ratio));
                exitProps = double(cea.get_exit_MolWt_gamma(Pc = pc, MR = MR, eps = expansion_ratio));
                chamberProps = double(cea.get_Chamber_MolWt_gamma(Pc = pc, MR = MR, eps = expansion_ratio));
                % lbm/lbm-mol -> kg/mol
                uniqueCEARuns.MW_t(u) = throatProps(1) * 1E-3;
                uniqueCEARuns.gamma_t(u) = throatProps(2);
                uniqueCEARuns.gamma_e(u) = exitProps(2);
                uniqueCEARuns.gamma_c(u) = chamberProps(2);

                uniqueCEARuns.pe(u) = pc / double(cea.get_PcOvPe(Pc = pc, MR = MR, eps = expansion_ratio));
                uniqueCEARuns.Isp_theo(u) = double(cea.get_Isp(Pc = pc, MR = MR, eps = expansion_ratio));

                % deg R to K
                temps = double(cea.get_Temperatures(Pc = pc, MR = MR, eps = expansion_ratio)) .* 5 ./ 9;
                uniqueCEARuns.T_c(u) = temps(1);
                uniqueCEARuns.T_t(u) = temps(2);
                uniqueCEARuns.T_e(u) = temps(3);
                % lbm/ft^3 -> kg/m^3
                densities = double(cea.get_Densities(Pc = pc, MR = MR, eps = expansion_ratio)) .* 16.018463;
                uniqueCEARuns.rho_t(u) = densities(2);

                % Transport properties: [cp, mu, k, Pr]
                throatTrans = double(cea.get_Throat_Transport(Pc = pc, MR = MR, eps = expansion_ratio));
                chamberTrans = double(cea.get_Chamber_Transport(Pc = pc, MR = MR, eps = expansion_ratio));
                % Btu/lbm-degR to J/kg-K
                uniqueCEARuns.cp_t(u) = throatTrans(1) * 1055.06 / (.453592 / 1.8);
                uniqueCEARuns.cp_c(u) = chamberTrans(1) * 1055.06 / (.453592 / 1.8);
                % milipoise to Pa-s
                uniqueCEARuns.mu_t(u) = throatTrans(2) * 1E-4;
                uniqueCEARuns.mu_c(u) = chamberTrans(2) * 1E-4;
                % millicalories / (cm-K-sec) to W/m-K
                uniqueCEARuns.k_t(u) = throatTrans(3) * .4184;
                uniqueCEARuns.k_c(u) = chamberTrans(3) * .4184;
                uniqueCEARuns.Pr_t(u) = throatTrans(4);
                uniqueCEARuns.Pr_c(u) = chamberTrans(4);

                % ft/s -> m/s
                uniqueCEARuns.cstar_theo(u) = double(cea.get_Cstar(Pc = pc, MR = MR)) * .3048;
            end

            % Spread back out to one value per engine
            for fieldName = fieldNames
                gasProps.(fieldName) = uniqueCEARuns.(fieldName)(mapToUnique);
            end
        end

        function geom = buildContour(in, dt, Vc)
        % Builds the nozzle and chamber contour of every engine at once
        %   r, z: N x S arrays, N engines, S stations from nozzle exit to
        %   injector
        % -----------------------------------------------------------------
        % Assumptions
        %   See Assumption 5 of the class header. Also:
        %   1) The nozzle's circular arc has radius .382 times the throat
        %       radius, the converging arc 1.5 times
        %   2) Chamber volume Vc is split into the converging sections and
        %       the straight section, by volume of revolution
        % -----------------------------------------------------------------
        % Sources
        %   1) Rao approximation bell nozzle, GT slides:
        %       https://www.seitzman.gatech.edu/classes/ae6450/nozzle_geometries.pdf
        % -----------------------------------------------------------------
            n = RocketEngine.numStations;
            Rt = dt ./ 2;
            thetaN = in.thetaN ./ 180 .* pi;

            % Sections are built about (0, 0) and shifted into place, from the
            % throat outward both ways. Order from nozzle exit to injector:
            % nozzle parabolic, nozzle circular, chamber circular, chamber
            % parabolic, chamber linear.

            %% Nozzle circular
            R1 = .382 .* Rt;
            xN = R1 .* sin(thetaN);
            L_nozzle_circular = xN;
            z_nozzle_circular = linspace(0, 1, n.nozzleCircular) .* L_nozzle_circular;
            r_nozzle_circular = (Rt + R1) - sqrt(R1.^2 - z_nozzle_circular.^2);
            % Flip so it runs from the parabolic section to the throat
            r_nozzle_circular = flip(r_nozzle_circular, 2);
            L_nozzle = in.fpct .* Rt ./ tan(thetaN) .* (sqrt(in.expansion_ratio) - 1 + 1.5 .* (1 ./ cos(thetaN) - 1));
            z_nozzle_circular = z_nozzle_circular + L_nozzle - L_nozzle_circular;

            %% Chamber circular
            R1p = Rt .* 1.5;
            Rc = sqrt(in.CR) .* Rt;
            % How far between throat radius and chamber radius the arc ends
            arcFraction = in.etaRc .* (Rc - Rt) ./ R1p;
            if any(arcFraction > 2)
                error("Converging arc is geometrically impossible for these CR and etaRc inputs, lower etaRc");
            end
            alpha = acos(1 - arcFraction);
            L_chamber_circular = R1p .* sin(alpha);
            z_chamber_circular = linspace(0, 1, n.chamberCircular) .* L_chamber_circular;
            r_chamber_circular = (Rt + R1p) - sqrt(R1p.^2 - z_chamber_circular.^2);
            z_chamber_circular = z_chamber_circular + z_nozzle_circular(:, end);

            %% Nozzle parabolic
            R_exit = sqrt(in.expansion_ratio) .* Rt;
            L_nozzle_parabolic = L_nozzle - L_nozzle_circular;
            z_nozzle_parabolic = linspace(0, 1, n.nozzleParabolic) .* L_nozzle_parabolic;
            % Parabola r = sqrt(b * (z - h)) + k, solved by hand from the
            % angle thetaN at the circular section and exit radius
            k = R_exit.^2 ./ (2 .* tan(thetaN) .* L_nozzle_parabolic + 2 .* R_exit);
            b = 2 .* tan(thetaN) .* k;
            h = -b ./ (4 .* tan(thetaN).^2);
            r_nozzle_parabolic = sqrt(b .* (z_nozzle_parabolic - h)) + k;
            % r isn't controlled to meet the circular section, so shift it to
            r_nozzle_parabolic = r_nozzle_parabolic + (r_nozzle_circular(:, 1) - r_nozzle_parabolic(:, 1));
            % Flip r so it runs from the exit plane in
            r_nozzle_parabolic = flip(r_nozzle_parabolic, 2);

            %% Chamber parabolic
            % Parabola r = a * z^2 + b * z that starts at the end of the
            % converging arc at its angle alpha, and ends tangent to the
            % chamber wall
            r1 = Rc - r_chamber_circular(:, end);
            b_p = tan(alpha);
            a_p = -b_p.^2 ./ (4 .* r1);
            L_chamber_parabolic = (-b_p + sqrt(max(b_p.^2 + 4 .* a_p .* r1, 0))) ./ (2 .* a_p);
            z_chamber_parabolic = linspace(0, 1, n.chamberParabolic) .* L_chamber_parabolic;
            r_chamber_parabolic = a_p .* z_chamber_parabolic.^2 + b_p .* z_chamber_parabolic;
            z_chamber_parabolic = z_chamber_parabolic + z_chamber_circular(:, end);
            r_chamber_parabolic = r_chamber_parabolic + r_chamber_circular(:, end);

            %% Chamber linear
            % What's left of the chamber volume after the converging sections
            Vcl = Vc - pi .* sum(r_chamber_circular(:, 1:end - 1).^2 .* diff(z_chamber_circular, 1, 2), 2);
            Vcl = Vcl - pi .* sum(r_chamber_parabolic(:, 1:end - 1).^2 .* diff(z_chamber_parabolic, 1, 2), 2);
            L_chamber_linear = Vcl ./ (pi .* Rc.^2);
            if any(L_chamber_linear < 0)
                warning("Chamber volume is smaller than the converging sections' volume for some engines (L* too small, or CR too large); their straight chamber length is negative");
            end
            r_chamber_linear = Rc .* ones(1, n.chamberLinear);
            z_chamber_linear = linspace(0, 1, n.chamberLinear) .* L_chamber_linear + max(z_chamber_parabolic, [], 2);

            %% Compile all sections
            geom.r = [r_nozzle_parabolic, r_nozzle_circular, r_chamber_circular, r_chamber_parabolic, r_chamber_linear];
            geom.z = [z_nozzle_parabolic, z_nozzle_circular, z_chamber_circular, z_chamber_parabolic, z_chamber_linear];
            geom.R1 = R1;
            geom.R1p = R1p;
            geom.Rc = Rc;
            geom.R_exit = R_exit;
            geom.thetaN = thetaN;
            geom.xN = xN;
            geom.alpha = alpha;
            geom.L_nozzle = L_nozzle;
            geom.L_nozzle_parabolic = L_nozzle_parabolic;
            geom.L_chamber_parabolic = L_chamber_parabolic;
            geom.L_chamber_linear = L_chamber_linear;
            % Station index of the throat (the last station of the nozzle
            % circular section. The first station of the chamber circular
            % section is a duplicate of it, at the same z)
            geom.iThroat = (n.nozzleParabolic + n.nozzleCircular) .* ones(size(Rt));
        end

        function vol = getShellVolume(geom, in)
        % Volume of the material of the engine (wall, channels and jacket)
        % from integrating the cross section over z                 vol in (m^3)
        % -----------------------------------------------------------------
        % Assumptions
        %   1) The jacket's outer radius is the inner radius plus 3 wall
        %       thicknesses, and the channels' cross section area (square,
        %       d_channel^2 each) is removed
        %   2) The earlier calculator summed this over stations without
        %       multiplying by station spacing, so its "volume" had units of
        %       area. This integrates (trapezoid rule). It also multiplied
        %       the channels' area by pi, which only suits round channels.
        % -----------------------------------------------------------------
            area = pi .* ((geom.r + 3 .* in.t_wall).^2 - geom.r.^2) - in.num_channels .* in.d_channel.^2;
            vol = sum((area(:, 1:end - 1) + area(:, 2:end)) ./ 2 .* diff(geom.z, 1, 2), 2);
        end

        function gasLine = getGasAlongContour(geom, gasProps)
        % Combustion gas Mach number, static temperature and transport
        % properties at every station
        % -----------------------------------------------------------------
        % Assumptions
        %   See Assumption 4 of the class header. Also:
        %   1) Static temperature uses the throat gas temperature to find
        %       the stagnation temperature, with the throat gamma, as the
        %       earlier calculator did: T = T_t * (1 + (gamma_t - 1)/2) /
        %       (1 + (gamma - 1)/2 * M^2), where gamma is the chamber's for
        %       subsonic stations and the throat's for the rest
        %   2) The gas velocity formula in the earlier calculator was
        %       missing the molar weight, so no velocity is stored. Mass
        %       flux, which is all the heat transfer calc uses, is mdot/A.
        % -----------------------------------------------------------------
            [~, numStations] = size(geom.r);
            Rt = geom.r(sub2ind(size(geom.r), (1:size(geom.r, 1))', geom.iThroat));
            areaRatio = geom.r.^2 ./ Rt.^2;

            % The last station of the throat region, counted from the nozzle
            % exit, is where subsonic flow stops
            nearThroat = abs(areaRatio - 1) <= RocketEngine.tol_throat;
            lastThroatStation = max(nearThroat .* (1:numStations), [], 2);
            subsonic = (1:numStations) > lastThroatStation;

            % Chamber values at subsonic stations, throat values everywhere
            % else
            full = ones(size(areaRatio));
            gamma_c = gasProps.gamma_c .* full;
            mu_c = gasProps.mu_c .* full;
            k_c = gasProps.k_c .* full;
            Pr_c = gasProps.Pr_c .* full;
            gamma_station = gasProps.gamma_t .* full;
            gamma_station(subsonic) = gamma_c(subsonic);
            gasLine.mu = gasProps.mu_t .* full;
            gasLine.mu(subsonic) = mu_c(subsonic);
            gasLine.k = gasProps.k_t .* full;
            gasLine.k(subsonic) = k_c(subsonic);
            gasLine.Pr = gasProps.Pr_t .* full;
            gasLine.Pr(subsonic) = Pr_c(subsonic);

            gasLine.AR = areaRatio;
            gasLine.subsonic = subsonic;
            gasLine.M = RocketEngine.GetMachFromAreaRatio(area_ratio = areaRatio, gamma = gamma_station, subsonic = subsonic);
            gasLine.T = gasProps.T_t .* (1 + (gasProps.gamma_t - 1) ./ 2) ./ (1 + (gamma_station - 1) ./ 2 .* gasLine.M.^2);
        end

        function coolant = getCoolantInlet(in, fuel)
        % Coolant state and flow in a channel at the channel inlet
        % -----------------------------------------------------------------
        % Assumptions
        %   See Assumptions 2 and 3 of the class header. Also:
        %   1) Channels are square, so area is d^2 and hydraulic diameter
        %       (4 * area / perimeter) is d
        % -----------------------------------------------------------------
            coolant.A_channel = in.d_channel.^2;
            coolant.dH = 4 .* coolant.A_channel ./ (4 .* in.d_channel);
            % Only fuel in the channels
            coolant.mdot_channel = in.mdot .* (1 ./ (1 + in.MR)) ./ in.num_channels;
            coolant.rho = fuel.getRho(in.T_coolant_i);
            coolant.k = fuel.getK(in.T_coolant_i, coolant.rho);
            coolant.mu = fuel.getMu(in.T_coolant_i, coolant.rho);
            coolant.cp = fuel.getCp(in.T_coolant_i);
            coolant.vel = coolant.mdot_channel ./ coolant.A_channel ./ coolant.rho;
            coolant.Re = coolant.dH .* coolant.vel .* coolant.rho ./ coolant.mu;
            coolant.Pr = coolant.cp .* coolant.mu ./ coolant.k;
        end

        function wall = getThroatEstimate(in, gasProps, perf, coolant)
        % Quick lumped estimate of heat flux and wall temperatures at the
        % throat. Cheaper than the station-wise solution
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Gas side: Dittus-Boelter type h_g from throat conditions,
        %       with the throat sonic flow velocity
        %   2) Coolant side: Dittus-Boelter for Re above Re_D_transition (n
        %       = .4, heated fluid), and the Colburn-type expression the
        %       earlier calculator used below it
        % -----------------------------------------------------------------
        % Sources
        %   1) https://www.sciencedirect.com/topics/engineering/dittus-boelter-correlation
        %   2) Equations 8-21 and 8-24 in Rocket Propulsion Elements by Sutton
        % -----------------------------------------------------------------
            Nu_coolant = .023 .* coolant.Re.^.8 .* coolant.Pr.^.4;
            h_l = Nu_coolant .* coolant.k ./ coolant.dH;
            laminar = coolant.Re <= RocketEngine.Re_D_transition;
            h_l_laminar = .023 .* coolant.cp .* coolant.mdot_channel ./ coolant.A_channel .* coolant.Re.^-.2 .* (coolant.mu .* coolant.cp ./ coolant.k).^(-2 ./ 3);
            h_l(laminar) = h_l_laminar(laminar);

            h_g = .023 .* (gasProps.rho_t .* perf.flow_sonic).^.8 ./ perf.dt.^.2 .* gasProps.Pr_t.^.4 .* gasProps.k_t ./ gasProps.mu_t.^.8;

            % Initial T_wg_avg guess
            T_wg_avg_guess = 0.5 * (gasProps.T_t + in.T_coolant_i);
            % Need k to get T_wg_throat, but need T_wg_throat to get k, so
            % iterate with T_wg_avg_guess until that temperature predicts
            for i = 1:RocketEngine.maxIter_T_wg
                k_wall = in.material.get_k(T_wg_avg_guess);
                TR_lumped = 1 ./ h_g + in.t_wall ./ k_wall + 1 ./ h_l ./ ((perf.dt + in.t_wall) ./ perf.dt);
                wall.q_throat = (gasProps.T_t - in.T_coolant_i) ./ TR_lumped;
                wall.T_wg_throat = gasProps.T_t - wall.q_throat ./ h_g;
                wall.T_wl_throat = in.T_coolant_i + wall.q_throat ./ h_l;
                if abs(T_wg_avg_guess - 0.5 * (wall.T_wg_throat + wall.T_wl_throat)) <= RocketEngine.tol_T_wg
                    break
                else
                    T_wg_avg_guess = 0.5 * (wall.T_wg_throat + wall.T_wl_throat);
                end
            end

            wall.h_g_throat = h_g;
            wall.h_l_throat = h_l;
        end

        function stations = solveWallTemperatures(in, material, gasProps, perf, geom, gasLine, coolant)
        % Heat flux and wall temperatures at every station of every engine,
        % by fixed-point iteration on the gas-side wall temperature, all
        % stations and engines at once
        %   stations: N x S arrays (and N x 1 for the time to steady state)
        % -----------------------------------------------------------------
        % Assumptions
        %   See Assumptions 3, 4 and 6 of the class header. Also:
        %   1) Hot wall k is a function of the mean wall temperature
        %   2) Gas-side h_g: Bartz, with chamber specific heat everywhere,
        %       throat gamma everywhere
        %   3) Coolant-side h_l: Gnielinski, using position along the
        %       channel measured from the nozzle exit, where coolant is
        %       taken to enter. (The old calculator subtracted the arc
        %       length from x again on every iteration of the temperature
        %       loop, which this doesn't)
        %   4) Time to steady state is the energy to heat the hot wall from
        %       T_amb to the gas temperature, over the steady heat flux
        %   5) Stations that haven't converged after maxIter_T_wg iterations
        %       keep their last iterate, flagged in converged.
        % -----------------------------------------------------------------
            [N, S] = size(geom.r);
            r = geom.r;
            t = in.t_wall;
            Tg = gasLine.T;
            A = gasLine.AR .* perf.At;

            % Local radius of curvature of the wall
            dzSeg = diff(geom.z, 1, 2);
            drSeg = diff(r, 1, 2);
            dz = [dzSeg(:, 1), dzSeg];
            dr = [drSeg(:, 1), drSeg];
            % Where the wall is too straight, set the radius of curvature to
            % the throat diameter, which effectively cancels Bartz's term
            R_curv = perf.dt .* ones(N, S);
            R_formula = RocketEngine.GetRadiusOfCurvature(dx = dr, dy = dz);
            % Sections meet at duplicate stations, so a segment can have no
            % length (just rounding noise), which has no curvature to speak of
            hasLength = sqrt(dr.^2 + dz.^2) > 1E-9 .* perf.dt;
            curved = abs(dr ./ dz) > RocketEngine.tol_R_curv & hasLength;
            R_curv(curved) = R_formula(curved);

            % Position along the channel, and total channel length
            segLength = sqrt(dzSeg.^2 + drSeg.^2);
            L_channel = sum(segLength, 2);
            x_channel = [zeros(N, 1), cumsum(segLength, 2)];

            % Coolant side film coefficient (doesn't depend on wall T)
            f = RocketEngine.GetHaalandFrictionFactor(epsilon = material.epsilon, D = coolant.dH, Re = coolant.Re);
            Nu_D = RocketEngine.GetGnielinskiNu_D(f = f, Re_D = coolant.Re, Pr = coolant.Pr, D = coolant.dH, L = L_channel, x = x_channel);
            h_l = RocketEngine.Geth(Nu_D = Nu_D, D = coolant.dH, k = coolant.k);
            TR_l = RocketEngine.GetThermalResistance(heat_transfer_type = "forced convection", h = h_l, r_0 = r, r_1 = r + t);

            T_wg_guess = .5 .* Tg;
            k_wall = material.get_k(T_wg_guess);

            active = true(N, S);
            T_wg = nan(N, S);
            T_wl = nan(N, S);
            q = nan(N, S);
            h_g_out = nan(N, S);
            TR_g_out = nan(N, S);
            for iter = 1:RocketEngine.maxIter_T_wg
                % Flow-wall convection
                sigma = RocketEngine.GetBartzCorrectionFactor(T_wg = T_wg_guess, T_c = Tg, gammay = gasProps.gamma_t, M = gasLine.M);
                h_g = RocketEngine.GetBartzh_g(D_t = perf.dt, mu = gasLine.mu, c_p = gasProps.cp_c, Pr = gasLine.Pr, p_c = in.pc, g = RocketEngine.g0, c_star = perf.cstar, R = R_curv, A_t = perf.At, A = A, sigma = sigma);
                TR_g = RocketEngine.GetThermalResistance(heat_transfer_type = "forced convection", h = h_g, r_0 = r, r_1 = r);

                % Hot wall conduction, first with k from the guess, then
                % again with k from the resulting mean wall temperature
                TR_w = RocketEngine.GetThermalResistance(heat_transfer_type = "conduction", k = k_wall, r_0 = r, r_1 = r, r_2 = r + t);
                q_calc = (Tg - in.T_coolant_i) ./ (TR_g + TR_w + TR_l);
                T_wg_calc = Tg - q_calc .* TR_g;
                T_wl_calc = in.T_coolant_i + q_calc .* TR_l;

                k_wall = material.get_k(.5 .* (T_wg_calc + T_wl_calc));
                TR_w = RocketEngine.GetThermalResistance(heat_transfer_type = "conduction", k = k_wall, r_0 = r, r_1 = r, r_2 = r + t);
                q_calc = (Tg - in.T_coolant_i) ./ (TR_g + TR_w + TR_l);
                T_wg_calc = Tg - q_calc .* TR_g;
                T_wl_calc = in.T_coolant_i + q_calc .* TR_l;

                % Check convergence, freeze the stations that have
                delta = T_wg_calc - T_wg_guess;
                justConverged = active & abs(delta) <= RocketEngine.tol_T_wg;
                T_wg(justConverged) = T_wg_calc(justConverged);
                T_wl(justConverged) = T_wl_calc(justConverged);
                q(justConverged) = q_calc(justConverged);
                h_g_out(justConverged) = h_g(justConverged);
                TR_g_out(justConverged) = TR_g(justConverged);
                active = active & ~justConverged;
                if ~any(active(:))
                    break;
                end
                T_wg_guess(active) = T_wg_guess(active) + .8 .* delta(active);
            end

            % Stations that never converged keep the last iterate
            if any(active(:))
                warning("Wall temperature didn't converge at %d of %d stations", nnz(active), numel(active));
                T_wg(active) = T_wg_calc(active);
                T_wl(active) = T_wl_calc(active);
                q(active) = q_calc(active);
                h_g_out(active) = h_g(active);
                TR_g_out(active) = TR_g(active);
            end

            stations.converged = ~active;
            stations.T_wg = T_wg;
            stations.T_wl = T_wl;
            stations.q_steady = q;
            stations.h_g = h_g_out;
            stations.h_l = h_l;
            % Start up: hot wall at T_amb, flux into it from the gas, and out
            % to the coolant
            stations.q_start_gas = (Tg - in.T_amb) ./ TR_g_out;
            stations.q_start_coolant = (in.T_amb - in.T_coolant_i) ./ TR_l;
            % Time to steady state
            delta_u = material.rho .* t .* material.cp .* (Tg - in.T_amb);
            time_to_steady_state = delta_u ./ q;
            stations.time_to_steady_state_min = min(time_to_steady_state, [], 2);
            stations.time_to_steady_state_max = max(time_to_steady_state, [], 2);
        end

        function merged = mergeStructs(a, b)
        % Adds the fields of struct b to struct a
            merged = a;
            fieldNames = fieldnames(b);
            for i = 1:numel(fieldNames)
                merged.(fieldNames{i}) = b.(fieldNames{i});
            end
        end

        function row = getRow(arrays, i)
        % Takes the i-th row of every field of a struct of arrays, to make
        % the struct for engine i
            fieldNames = fieldnames(arrays);
            for j = 1:numel(fieldNames)
                row.(fieldNames{j}) = arrays.(fieldNames{j})(i, :);
            end
        end

    end

end

% -------------------------------------------------------------------------
% Method template, updated from the old FunctionTemplate. What was wrong
% with the old one, for a method of a class:
%   1) The first input of a non-static method is the object, so the
%       "args" name-value struct has to come second, or the arguments
%       block errors. Static methods (no object) are fine as they were.
%   2) eval() to make variables out of the args struct hides the variables
%       from MATLAB's code analyzer, is slow, and can't tell an argument
%       that was left out from one that was passed as [], so now the
%       required-argument check is RocketEngine.requireArgs, and the
%       method reads args.<name> into local variables with plain
%       assignments.
%
%         function result = MethodTemplate(obj, args)
%         % <Method Purpose>
%         %                                             <Output> in (<Units>)
%         % -----------------------------------------------------------------
%         % Arguments
%         %   <Symbol> = <Explanation> (<Units>)
%         % -----------------------------------------------------------------
%         % Dependencies
%         %   #) <Dependency Filepath>
%         % -----------------------------------------------------------------
%         % Assumptions
%         %   #) <Assumption>, or: See Assumption # of the class header
%         % -----------------------------------------------------------------
%         % Sources
%         %   #) <Source>
%         % -----------------------------------------------------------------
%         % MATLAB Version <Oldest Version>, also compatible with:
%         %   - <Later Version>
%         % -----------------------------------------------------------------
%             % Allows arguments to be assigned in the method call as in:
%             %   obj.MethodTemplate(<arg_name> = <arg_val>, ...)
%             arguments
%                 obj
%                 args.arg_1 = [];
%             end
%             % Names of the arguments that are not optional
%             RocketEngine.requireArgs(args, "arg_1");
%             arg_1 = args.arg_1;
%
%             % Unit Conversions
%
%             % Intermediate Calculations
%
%             % Final Calculations
%
%             % Display Results and/or Plotting
%
%             return;
%         end