%%
GSI_model1 = @coeff_DRIA;
GSI_model2 = @coeff_CLL;

h = 300e3;
lat = 0;
lon = 0;
dayOfYear = 1;
UTseconds = 0;
f107Average = 140;
f107Daily = 140;
magneticIndex = ones([1,7])*15;

param_eq = struct;
param_eq = environment(param_eq, h, lat, lon, dayOfYear, UTseconds, f107Average, f107Daily, magneticIndex, 1);

param_eq.Tw = 300;
delta = deg2rad(0:1:90); % Angle from Normal

param_eq.alpha = 1;
param_eq.gamma = cos(delta);
param_eq.ell = sin(delta);
[cp1, ctau1, cd1, cl1] = GSI_model1(param_eq, delta);

% CLL
param_eq.alphaN = 0.5;
param_eq.sigmaT = 0.5;
[cp2, ctau2, cd2, cl2] = GSI_model2(param_eq, delta);

figure
figuresize(18, 9, 'cm')
hold on
p_1 = plot(delta*(180/pi),cd1);
set(p_1, 'Color', 'k', 'LineStyle','--', 'LineWidth',1.25)
p_2 = plot(delta*(180/pi),cd2);
set(p_2, 'Color', 'k', 'LineStyle','-','LineWidth',1.25)
grid on
ylabel('Drag Coefficient, $C_D$')

yyaxis right
set(gca, 'YColor',[0.5 0.5 0.5])
hold on
p_3 = plot(delta*(180/pi),cl1);
set(p_3, 'Color', [0.5 0.5 0.5], 'LineStyle','--','LineWidth',1.25)
p_4 = plot(delta*(180/pi),cl2);
set(p_4, 'Color', [0.5 0.5 0.5], 'LineStyle','-','LineWidth',1.25)
ylabel('Lift Coefficient, $C_L$')

xlabel('Incidence Angle [deg]')
legend('$C_D$ DRIA ($\alpha = 1$)', '$C_D$ CLL ($\alpha_N = \alpha_T = 0.5$)', '$C_L$ DRIA ($\alpha = 1$)', '$C_L$ CLL ($\alpha_N = \alpha_T = 0.5$)')

set(gcf,'color','w');