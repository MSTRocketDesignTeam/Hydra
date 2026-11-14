% RocketEngine debug script

clc
clear


var_range = [1]';
T_amb = 21 + 273.15;

% Sweep these
p_c_array = [128] * var_range * 6894.76; % psi to Pa
MR_array = [1.6] * var_range;
mdot_array = [1.2 .* .45359237] * var_range; % kg/s
d_channel_array = [1.2] * 1E-3 * var_range; % m
num_channels_array = [48] * var_range;
T_coolant_i_array = [T_amb] * var_range;
t_wall_array = [0.8] * var_range *1E-3; % m
CR_array = [6] * var_range;

% one input per
name = "Hydra";
fuel = Propellant("ETHANOL_L");
oxidizer = Propellant("N2O_L");
material = Material("AlSi10Mg");
alt = 342;
Lstar = 0.6;
thetaN = 20;
fpct = 1;
etaRc = 0.2;
DF_stress = 2;
DF_temp = 1.5;
DF_p_coolant = 1.2;

% make RocketEngine object
engines = RocketEngine(name = name, ...
    fuel = fuel, ...
    oxidizer = oxidizer, ...
    material = material, ...
    alt = alt, ...
    T_amb = T_amb, ...
    p_c = p_c_array, ...
    MR = MR_array, ...
    mdot = mdot_array, ...
    d_channel = d_channel_array, ...
    num_channels = num_channels_array, ...
    T_coolant_i = T_coolant_i_array, ...
    t_wall = t_wall_array, ...
    CR = CR_array, ...
    Lstar = Lstar, ...
    thetaN = thetaN, ...
    fpct = fpct, ...
    etaRc = etaRc, ...
    DF_stress = DF_stress, ...
    DF_temp = DF_temp, ...
    DF_p_coolant = DF_p_coolant)

engine = engines(1);
thetaArray = linspace(0, 2 * pi, 20);
zArray = engine.geom.z .* 1000;
rArray = engine.geom.r .* 1000;
[thetaArray, zArray] = meshgrid(thetaArray, zArray);
rArray = repmat(rArray(:), 1, size(thetaArray, 2));
xArray = rArray .* cos(thetaArray);
yArray = rArray .* sin(thetaArray);
close all
fig = figure();
ax = axes(fig);
surf(ax, xArray, yArray, zArray);
axis equal;
view(3);