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
%       vectorizable, so it is called only for each unique (pc, MR, CR,
%       p_amb, thetaN, fpct, etaRc), in a few small calls each (see
%       runCEA).
%   10) The nozzle is ideally expanded at the design altitude: the exit
%       pressure equals the ambient pressure p_amb. CEA is given the
%       pressure ratio pc / p_amb and returns the expansion ratio, so
%       expansion_ratio is an output, not an input. pc is the injector face
%       pressure (CEA's Pinj), the small pressure drop to the end of the
%       chamber is kept.
%   11) The exhaust is frozen from the throat on (CEA frozen = 1 and
%       frozenAtThroat = 1). Chamber and throat state is equilibrium.
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
        %
        % CEA is asked for the gas state at the area ratio of every station.
        % The distinct area ratios it must report are
        %   supersonic: nozzleParabolic + nozzleCircular - 2 (the throat and
        %       the duplicate station where the two sections meet; the exit
        %       is one of them)
        %   subsonic: chamberCircular + chamberParabolic - 3 (the throat,
        %       the duplicate where the sections meet, and the chamber end)
        % Limits that matter (read from RocketCEA's Fortran): a CEA call
        % takes at most 8 supar and, with a finite area combustor, 12 subar,
        % and a frozen-at-throat call silently drops its frozen results if it
        % has more than 3 extra points. runCEA splits the lists across
        % calls (maxSuparPerCall, maxSubarPerCall), so there is no hard cap
        % on these numbers. Every 3 more supersonic or 4 more subsonic
        % stations is one more CEA call per (pc, MR) case.
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

        % CEA call splitting, see numStations
        maxSuparPerCall = 3;
        maxSubarPerCall = 4;
        % Station area ratios closer than this (relative) are one CEA point.
        % Half of it is also how closely CEA's printed area ratios must match
        % the requested ones
        tol_arSchedule = 1E-3;
        % Transport properties (cp, k, Pr) at the chamber, throat and
        % subsonic stations: false = equilibrium reactions (what the earlier
        % calculator used), true = frozen reactions. Supersonic stations
        % are always frozen, the exhaust is frozen from the throat on.
        useFrozenTransport = false;

        % Names of the inputs that can be swept, in the order of the dimensions
        % of the engine array
        sweepNames = ["pc", "MR", "mdot", "d_channel", ...
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
        expansion_ratio = []; % Nozzle exit area to throat area, computed by solve for exit pressure = ambient pressure
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
        %   pc, MR, mdot, d_channel, num_channels,
        %   T_coolant_i, t_wall, CR, Lstar = see properties.
        %       (expansion_ratio is not an input, solve computes it)
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
            % Also finds the expansion ratio (added to inputs) and the gas
            % state at every contour station
            [gasProps, gasLine, inputs] = RocketEngine.runCEA(inputs, fuel, oxidizer);

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
                objs(i).expansion_ratio = inputs.expansion_ratio(i);
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

        function cases = parseCEAOutput(outputText)
        % Reads the text from RocketCEA's get_full_cea_output (short output)
        % into a struct array, one element per case (a pc, O/F pair)
        %                                  cases(i).eq, equilibrium page
        %                                  cases(i).fr, frozen page (empty if
        %                                  the call wasn't frozen)
        % -----------------------------------------------------------------
        % Each page is a struct of row vectors with one value per column of
        % the page (injector, chamber end, throat, then the other points, in
        % CEA's units): PinjOvP, P, T, rho, MW, gamma, mach, mu, cpThermo,
        % cpEq, kEq, PrEq, cpFr, kFr, PrFr, and AeAt, cstar, CF, Ivac, Isp
        % (CEA doesn't print these for the injector, so their first value is
        % NaN), and OF, Pinj. Rows a page doesn't have are left out.
        % -----------------------------------------------------------------
        % Assumptions
        %   1) Every page fits in one block of columns (CEA breaks wider
        %       output into several blocks). The calls from runCEA are
        %       kept narrow enough, and this errors if a row's length
        %       isn't the page's number of columns
        % -----------------------------------------------------------------
            lines = regexp(char(outputText), '\r?\n', 'split');
            pageStarts = find(contains(lines, "THEORETICAL ROCKET PERFORMANCE ASSUMING"));
            if isempty(pageStarts)
                shown = char(outputText);
                error("No CEA results found in the CEA output. It starts: %s", shown(1:min(300, end)));
            end
            pageEnds = [pageStarts(2:end) - 1, numel(lines)];
            cases = struct("eq", {}, "fr", {});
            for p = 1:numel(pageStarts)
                page = RocketEngine.parseCEAPage(lines(pageStarts(p):pageEnds(p)));
                if page.isEquilibrium
                    cases(end + 1).eq = page; %#ok<AGROW>
                    cases(end).fr = [];
                else
                    if isempty(cases) || ~isempty(cases(end).fr)
                        error("CEA output has a frozen page that doesn't follow an equilibrium page");
                    end
                    cases(end).fr = page;
                end
            end
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

            inputNames = ["alt", "T_amb", "p_amb", "pc", "MR", ...
                "mdot", "d_channel", "num_channels", "T_coolant_i", "t_wall", "CR", "Lstar", "thetaN", "fpct", "etaRc", ...
                "DF_stress", "DF_temp", "DF_p_coolant"];
            for inputName = inputNames
                inputs.(inputName) = reshape([objs.(inputName)], [], 1);
            end
        end

        function [gasProps, gasLine, inputs] = runCEA(inputs, fuel, oxidizer)
        % Gets combustion gas state from CEA for every engine, in three
        % stages, calling CEA only for each unique (pc, MR, CR, p_amb,
        % thetaN, fpct, etaRc):
        %   1) Small calls with the pressure ratio pc / p_amb, which return
        %       the expansion ratio for exit pressure = ambient pressure
        %   2) Area ratio of every contour station, from a unit throat
        %       contour (area ratio doesn't depend on the throat size, so
        %       no resizing is needed after the throat is sized)
        %   3) Full calls at those area ratios, split into several calls
        %       each (see numStations)
        %                  gasProps, one value per engine
        %                  gasLine, N x S arrays of state at every station
        %                  inputs, the same with expansion_ratio added
        % -----------------------------------------------------------------
        % Assumptions
        %   1) CEA's inputs and outputs are English units, converted here
        %   2) See Assumptions 10 and 11 of the class header. The frozen
        %       exit gamma replaces the equilibrium one of the earlier
        %       version in the thrust coefficient.
        %   3) CEA's tables are matched to stations by area ratio, the
        %       nearest one (stations closer than tol_arSchedule to one
        %       another share a CEA point)
        % -----------------------------------------------------------------
            % Import RocketCEA library
            py.importlib.import_module("rocketcea");

            numEngines = numel(inputs.pc);
            psiPerPa = 0.000145038;
            atmToPa = 101325;

            CEAkeys = [inputs.pc, inputs.MR, inputs.CR, inputs.p_amb, ...
                inputs.thetaN, inputs.fpct, inputs.etaRc];
            [uniqueCEAKeys, firstIdx, mapToUnique] = unique(CEAkeys, "rows", "first");
            numUnique = size(uniqueCEAKeys, 1);

            % One CEA object per contraction ratio, since CR is fixed when
            % the object is made
            ceaByCR = containers.Map("KeyType", "double", "ValueType", "any");
            for CR = unique(uniqueCEAKeys(:, 3))'
                ceaByCR(CR) = py.rocketcea.cea_obj.CEA_Obj(fuelName = fuel.name, oxName = oxidizer.name, fac_CR = CR);
            end

            %% Stage 1: expansion ratio of each case
            % One call per (pc, CR, p_amb), with all of its MR's in a list
            epsUnique = nan(numUnique, 1);
            [uniqueGroups, ~, groupOf] = unique(uniqueCEAKeys(:, [1 3 4]), "rows");
            for g = 1:size(uniqueGroups, 1)
                members = find(groupOf == g);
                pc = uniqueGroups(g, 1);
                CR = uniqueGroups(g, 2);
                p_amb = uniqueGroups(g, 3);
                MRs = uniqueCEAKeys(members, 2);
                cases = RocketEngine.callCEA(ceaByCR(CR), pc * psiPerPa, MRs, [], [], pc / p_amb, 1, 0);
                if numel(cases) ~= numel(members)
                    error("CEA returned %d cases for %d O/F's (pc = %g Pa, CR = %g)", numel(cases), numel(members), pc, CR);
                end
                for m = 1:numel(members)
                    frozenBlock = cases(m).fr;
                    if abs(cases(m).eq.OF - MRs(m)) > 1E-4 .* MRs(m) || isempty(frozenBlock) || numel(frozenBlock.T) ~= 4
                        error("CEA's expansion ratio call for pc = %g Pa, O/F = %g, p_amb = %g Pa wasn't as expected. Is p_amb above the chamber end pressure?", pc, MRs(m), p_amb);
                    end
                    % The last column is the exit, at pc / p_amb
                    epsUnique(members(m)) = frozenBlock.AeAt(end);
                end
            end
            inputs.expansion_ratio = epsUnique(mapToUnique);

            %% Stage 2: area ratio at every contour station
            % Contour for a throat of radius 1, only to read area ratios
            unitGeom = RocketEngine.buildContour(inputs, 2 .* ones(numEngines, 1), 1E3 .* ones(numEngines, 1));
            [~, numStationsTotal] = size(unitGeom.r);
            rThroat = unitGeom.r(sub2ind(size(unitGeom.r), (1:numEngines)', unitGeom.iThroat));
            areaRatio = (unitGeom.r ./ rThroat).^2;

            %% Stage 3: full CEA calls
            if RocketEngine.useFrozenTransport
                subsonicTransport = "fr";
            else
                subsonicTransport = "eq";
            end
            gasPropFields = ["MW_t", "gamma_t", "gamma_e", "gamma_c", "pe", ...
                "Isp_theo", "T_c", "T_t", "T_e", "rho_t", "cp_t", "mu_t", ...
                "k_t", "Pr_t", "cp_c", "mu_c", "k_c", "Pr_c", "cstar_theo"];
            stateFields = ["M", "T", "P", "gamma", "mu", "k", "Pr", "cp", "rho"];
            gasPropsUnique = struct();
            for fieldName = gasPropFields
                gasPropsUnique.(fieldName) = zeros(numUnique, 1);
            end
            gasLineUnique = struct();
            for fieldName = stateFields
                gasLineUnique.(fieldName) = zeros(numUnique, numStationsTotal);
            end

            for u = 1:numUnique
                pc = uniqueCEAKeys(u, 1);
                MR = uniqueCEAKeys(u, 2);
                CR = uniqueCEAKeys(u, 3);
                cea = ceaByCR(CR);
                e = firstIdx(u);
                iThroat = unitGeom.iThroat(e);
                stationAR = areaRatio(e, :);

                % Distinct area ratios to ask CEA for
                supList = RocketEngine.getAreaRatioSchedule(stationAR(1:iThroat), 1, Inf);
                subList = RocketEngine.getAreaRatioSchedule(stationAR(iThroat + 1:end), 1, CR);
                if isempty(supList)
                    error("No supersonic stations found for engine %d", e);
                end

                % Supersonic stations, frozen at the throat. The first call
                % also gives the chamber and throat (equilibrium)
                supStates = [];
                supChunks = RocketEngine.splitIntoChunks(supList, RocketEngine.maxSuparPerCall);
                for c = 1:numel(supChunks)
                    cases = RocketEngine.callCEA(cea, pc * psiPerPa, MR, supChunks{c}, [], [], 1, 1);
                    if numel(cases) ~= 1 || isempty(cases(1).fr)
                        error("CEA didn't return one frozen case for pc = %g Pa, O/F = %g. CEA drops the frozen results when a call has too many points, see numStations.", pc, MR);
                    end
                    frozenBlock = cases(1).fr;
                    cols = RocketEngine.matchColumns(frozenBlock, 4, supChunks{c});
                    supStates = RocketEngine.catStates(supStates, RocketEngine.getCEAState(frozenBlock, cols, "fr", atmToPa));
                    if c == 1
                        % Chamber end is column 2 and throat column 3
                        chamberState = RocketEngine.getCEAState(cases(1).eq, 2, subsonicTransport, atmToPa);
                        throatState = RocketEngine.getCEAState(cases(1).eq, 3, subsonicTransport, atmToPa);
                        cstar = cases(1).eq.cstar(2) .* 0.3048;
                    end
                end

                % Subsonic stations, equilibrium
                subStates = [];
                subChunks = RocketEngine.splitIntoChunks(subList, RocketEngine.maxSubarPerCall);
                for c = 1:numel(subChunks)
                    cases = RocketEngine.callCEA(cea, pc * psiPerPa, MR, [], subChunks{c}, [], 0, 1);
                    if numel(cases) ~= 1
                        error("CEA returned %d cases for the subsonic area ratios at pc = %g Pa, O/F = %g", numel(cases), pc, MR);
                    end
                    cols = RocketEngine.matchColumns(cases(1).eq, 4, subChunks{c});
                    subStates = RocketEngine.catStates(subStates, RocketEngine.getCEAState(cases(1).eq, cols, subsonicTransport, atmToPa));
                end

                % Tables of state against area ratio, throat first
                nozzleTable = RocketEngine.catStates(throatState, supStates);
                nozzleTable.ar = [1; supList(:)];
                chamberTable = RocketEngine.catStates(throatState, subStates);
                chamberTable = RocketEngine.catStates(chamberTable, chamberState);
                chamberTable.ar = [1; subList(:); CR];

                % Spread onto the stations
                nozzleIdx = 1:iThroat;
                chamberIdx = iThroat + 1:numStationsTotal;
                nozzleLookup = RocketEngine.nearestIndex(nozzleTable.ar, stationAR(nozzleIdx));
                chamberLookup = RocketEngine.nearestIndex(chamberTable.ar, stationAR(chamberIdx));
                for fieldName = stateFields
                    gasLineUnique.(fieldName)(u, nozzleIdx) = nozzleTable.(fieldName)(nozzleLookup)';
                    gasLineUnique.(fieldName)(u, chamberIdx) = chamberTable.(fieldName)(chamberLookup)';
                end

                % Single values
                gasPropsUnique.MW_t(u) = throatState.MW;
                gasPropsUnique.gamma_t(u) = throatState.gamma;
                gasPropsUnique.gamma_e(u) = supStates.gamma(end);
                gasPropsUnique.gamma_c(u) = chamberState.gamma;
                gasPropsUnique.pe(u) = supStates.P(end);
                gasPropsUnique.Isp_theo(u) = supStates.Ivac(end);
                gasPropsUnique.T_c(u) = chamberState.T;
                gasPropsUnique.T_t(u) = throatState.T;
                gasPropsUnique.T_e(u) = supStates.T(end);
                gasPropsUnique.rho_t(u) = throatState.rho;
                gasPropsUnique.cp_t(u) = throatState.cp;
                gasPropsUnique.mu_t(u) = throatState.mu;
                gasPropsUnique.k_t(u) = throatState.k;
                gasPropsUnique.Pr_t(u) = throatState.Pr;
                gasPropsUnique.cp_c(u) = chamberState.cp;
                gasPropsUnique.mu_c(u) = chamberState.mu;
                gasPropsUnique.k_c(u) = chamberState.k;
                gasPropsUnique.Pr_c(u) = chamberState.Pr;
                gasPropsUnique.cstar_theo(u) = cstar;
            end

            % Spread back out to one value per engine
            for fieldName = gasPropFields
                gasProps.(fieldName) = gasPropsUnique.(fieldName)(mapToUnique);
            end
            for fieldName = stateFields
                gasLine.(fieldName) = gasLineUnique.(fieldName)(mapToUnique, :);
            end
            gasLine.AR = areaRatio;
            gasLine.subsonic = (1:numStationsTotal) > unitGeom.iThroat;
        end

        function cases = callCEA(cea, pcPsi, MR, supar, subar, PcOvPe, frozen, showTransport)
        % One get_full_cea_output call, parsed. Empty supar, subar or PcOvPe
        % are left out of the CEA input (None).
        %   frozen, 1 = frozen after the throat (frozen = 1 and
        %       frozenAtThroat = 1), 0 = equilibrium throughout
        % -----------------------------------------------------------------
        % Assumptions
        %   1) cea is a variable, not a chained call: chaining the method
        %       onto ceaByCR(CR) passes the name = value arguments wrong
        % -----------------------------------------------------------------
            optional = {supar, subar, PcOvPe};
            for j = 1:numel(optional)
                if isempty(optional{j})
                    optional{j} = py.None;
                else
                    optional{j} = double(optional{j}(:)');
                end
            end
            outputText = cea.get_full_cea_output(Pc = pcPsi, MR = double(MR(:)'), ...
                eps = optional{1}, subar = optional{2}, PcOvPe = optional{3}, ...
                frozen = frozen, frozenAtThroat = frozen, short_output = 1, ...
                show_transport = showTransport);
            cases = RocketEngine.parseCEAOutput(char(outputText));
        end

        function schedule = getAreaRatioSchedule(areaRatios, lowerLimit, upperLimit)
        % Distinct area ratios, ascending, strictly between the limits
        % (those within tol_arSchedule of a limit are left out). Area
        % ratios within tol_arSchedule of each other are merged, keeping the
        % larger one (so the nozzle exit is kept)
            tol = RocketEngine.tol_arSchedule;
            candidates = sort(areaRatios(:), "descend");
            candidates = candidates(candidates > lowerLimit .* (1 + tol) & candidates < upperLimit .* (1 - tol));
            schedule = zeros(0, 1);
            for k = 1:numel(candidates)
                if isempty(schedule) || candidates(k) < schedule(end) .* (1 - tol)
                    schedule(end + 1, 1) = candidates(k); %#ok<AGROW>
                end
            end
            schedule = flip(schedule);
        end

        function chunks = splitIntoChunks(values, maxPerChunk)
        % Splits a list into a cell array of consecutive pieces of at most
        % maxPerChunk
            values = values(:);
            numChunks = ceil(numel(values) ./ maxPerChunk);
            chunks = cell(1, numChunks);
            for c = 1:numChunks
                chunks{c} = values((c - 1) .* maxPerChunk + 1:min(c .* maxPerChunk, numel(values)));
            end
        end

        function cols = matchColumns(block, firstColumn, requested)
        % Columns of a CEA block (from firstColumn on) whose area ratio is
        % each requested one, in the order requested
            candidates = firstColumn:numel(block.T);
            candidateAR = block.AeAt(candidates);
            requested = requested(:)';
            [relErr, k] = min(abs(candidateAR(:) - requested) ./ requested, [], 1);
            if numel(candidates) < numel(requested) || any(relErr > RocketEngine.tol_arSchedule ./ 2) || numel(unique(k)) < numel(k)
                error("CEA's area ratios (%s) don't match the requested ones (%s)", ...
                    num2str(candidateAR(:)', "%g "), num2str(requested, "%g "));
            end
            cols = candidates(k);
        end

        function idx = nearestIndex(tableValues, queryValues)
        % For each query value, the index of the nearest table value
            [~, idx] = min(abs(queryValues(:) - tableValues(:)'), [], 2);
        end

        function states = catStates(a, b)
        % Stacks two state structs (see getCEAState) into one
            if isempty(a)
                states = b;
                return;
            elseif isempty(b)
                states = a;
                return;
            end
            states = a;
            for fieldName = string(fieldnames(b))'
                states.(fieldName) = [a.(fieldName); b.(fieldName)];
            end
        end

        function state = getCEAState(block, cols, transport, atmToPa)
        % Gas state at the columns of a parsed CEA block, in SI
        %   transport, "eq" or "fr": which cp, k and Pr rows to use
            if transport == "eq"
                cp = block.cpEq;
                k = block.kEq;
                Pr = block.PrEq;
            else
                cp = block.cpFr;
                k = block.kFr;
                Pr = block.PrFr;
            end
            state.M = block.mach(cols)';
            state.T = block.T(cols)';
            % atm to Pa
            state.P = block.P(cols)' .* atmToPa;
            % g/cc to kg/m^3
            state.rho = block.rho(cols)' .* 1000;
            % g/mol to kg/mol
            state.MW = block.MW(cols)' .* 1E-3;
            state.gamma = block.gamma(cols)';
            % millipoise to Pa-s
            state.mu = block.mu(cols)' .* 1E-4;
            % cal/g-K to J/kg-K
            state.cp = cp(cols)' .* 4184;
            % millicalories / (cm-K-sec) to W/m-K
            state.k = k(cols)' .* .4184;
            state.Pr = Pr(cols)';
            % s (CEA's LB-SEC/LB)
            state.Ivac = block.Ivac(cols)';
        end

        function page = parseCEAPage(pageLines)
        % Reads one page of CEA's output, see parseCEAOutput
            page = struct();
            page.isEquilibrium = contains(pageLines{1}, "ASSUMING EQUILIBRIUM");
            rowKeys = ["PinjOvP", "P", "T", "rho", "MW", "gamma", "mach", "mu", ...
                "cp", "k", "Pr", "AeAt", "cstar", "CF", "Ivac", "Isp"];
            rowLabels = ["Pinj/P", "P, ATM", "T, K", "RHO, G/CC", "MW, MOL WT", ...
                "GAMMAs", "MACH NUMBER", "VISC,MILLIPOISE", "Cp, CAL/(G)(K)", ...
                "CONDUCTIVITY", "PRANDTL NUMBER", "Ae/At", "CSTAR, FT/SEC", "CF", ...
                "Ivac,LB-SEC/LB", "Isp, LB-SEC/LB"];
            % The transport rows come in sections
            section = "thermo";
            for i = 1:numel(pageLines)
                text = strtrim(pageLines{i});
                if startsWith(text, "Pinj =")
                    page.Pinj = str2double(regexp(text, '[-+]?\d+\.?\d*', 'match', 'once'));
                elseif startsWith(text, "O/F=")
                    page.OF = str2double(regexp(text, '[-+]?\d+\.?\d*', 'match', 'once'));
                elseif startsWith(text, "WITH EQUILIBRIUM REACTIONS")
                    section = "eq";
                elseif startsWith(text, "WITH FROZEN REACTIONS")
                    section = "fr";
                elseif startsWith(text, "PERFORMANCE PARAMETERS")
                    section = "perf";
                elseif startsWith(text, "MOLE FRACTIONS")
                    break;
                else
                    for j = 1:numel(rowLabels)
                        if startsWith(text, rowLabels(j))
                            values = RocketEngine.getCEAValues(extractAfter(text, strlength(rowLabels(j))));
                            % Lines like "CONDUCTIVITY IN UNITS OF" have no values
                            if isempty(values)
                                break;
                            end
                            key = rowKeys(j);
                            if any(key == ["cp", "k", "Pr"])
                                if key == "cp" && section == "thermo"
                                    key = "cpThermo";
                                else
                                    key = key + upper(extractBefore(section, 2)) + extractAfter(section, 1);
                                end
                            end
                            page.(key) = values;
                            break;
                        end
                    end
                end
            end

            if ~isfield(page, "T") || ~isfield(page, "P") || ~isfield(page, "AeAt")
                error("A page of the CEA output is missing its temperature, pressure or area ratio row");
            end
            numColumns = numel(page.T);
            % Performance rows skip the injector column, so pad them
            for key = ["AeAt", "cstar", "CF", "Ivac", "Isp"]
                if isfield(page, key)
                    page.(key) = [NaN, page.(key)];
                end
            end
            for key = string(fieldnames(page))'
                if isnumeric(page.(key)) && numel(page.(key)) > 1 && numel(page.(key)) ~= numColumns
                    error("CEA output row %s has %d values on a page of %d columns. The page may be split into blocks of columns.", key, numel(page.(key)), numColumns);
                end
            end
        end

        function values = getCEAValues(text)
        % The numbers in a line of CEA output. CEA writes small ones without
        % the E, like 1.3106-3
            tokens = regexp(text, '[-+]?\d+\.\d+(?:[-+]\d+)?', 'match');
            values = zeros(1, numel(tokens));
            for i = 1:numel(tokens)
                parts = regexp(tokens{i}, '^([-+]?\d+\.\d+)([-+]\d+)$', 'tokens', 'once');
                if isempty(parts)
                    values(i) = str2double(tokens{i});
                else
                    values(i) = str2double(parts{1}) .* 10 .^ str2double(parts{2});
                end
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
        % NOT CALLED ANYMORE: runCEA gets the gas state at every station from
        % CEA itself. Kept as the ideal gas version, as a cross check.
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