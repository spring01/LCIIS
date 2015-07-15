% cart = [...
%  8                 -2.54062035   -0.22895125    0.00000000
%  1                 -1.58062035   -0.22895125    0.00000000
%  1                 -2.86107494    0.67598458    0.00000000];

cart = [...
    8 0.0 0.0 0.0
    1 0.0 1.0 0.0
    1 0.0 0.0 1.0];

% cart = [...
%  14                -1.55096010   -0.68685376    0.00000000
%  1                 -1.06097692   -2.07278899    0.00000000
%  1                 -1.06095162    0.00610450    1.20025020
%  1                 -1.06095162    0.00610450   -1.20025020
%  1                 -5.55096010   -0.68680447    0.00000000];

% cart = [...
%     14                  0.000000    0.000000    0.140556
%     1                   0.000000    1.385929    0.630556
%     1                   1.200250   -0.692965    0.630556
%     1                  -1.200250   -0.692965    0.630556
%     1                   0.000000    0.000000   -3.859444];

% cart = [...
%  24                0.00000000    0.00000000    0.40000000
%  6                 0.00000000    0.00000000   -1.60000000];


% cart = [...
%     48     0.000000     0.000000     0.000000
%     7      0.000000     0.000000    -2.260001
%     7     -0.685444     0.000000    -4.348035
%     6      0.676053     0.000000    -4.385069
%     6      1.085240     0.000000    -3.091231
%     6     -1.044752     0.000000    -3.060220
%     1      1.231530     0.000000    -5.300759
%     1      2.088641     0.000000    -2.711077
%     1     -2.068750     0.000000    -2.726515
%     1     -1.313170     0.000000    -5.174718];

% cart = [...
% 6                  -1.171138   -0.148546    0.000004
% 1                  -1.714145    0.219764   -0.880470
% 1                  -1.713691    0.218958    0.881075
% 1                  -1.155613   -1.240682   -0.000601
% 6                   0.234078    0.399597    0.000022
% 1                   0.301980    1.511628   -0.000146
% 8                   1.237979   -0.276996   -0.000002];

import SCF.*;

mol = Molecule(cart);
basisSet = '6-31g*';
dft = 'b3lyp';
diisType = 'C20';

matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
% matpsi.Settings_SetMaxNumCPUCores(2);
matpsi.SCF_SetSCFType('uhf');
matpsi.JK_Initialize('directjk');

% scf = RHF(RHF.MatPsi2Interface(matpsi));
scf = UKS(RHF.MatPsi2Interface(matpsi), dft);
[guessDensity, guessOrbital] = scf.CoreGuess();

[ener1, energySet1, iter1] = scf.SCF(guessOrbital, diisType);

% info.chargeMult = [0 1];
% info.cartesian = cart;
% info.method = 'b3lyp';
% info.basisSet = '6-31g';
% scf2 = G09RSCF(info);
% [guessDensity, guessOrbital] = scf2.CoreGuess();
% [ener1, energySet1, iter1] = scf2.SCF(guessOrbital, diisType);

fprintf('%0.8f  %d \n',ener1, iter1);

figure();
hold();
plot(log10(abs(energySet1 - ener1)), 'r');
scatter(1:length(energySet1), log10(abs(energySet1 - ener1)), 72, 'square', 'r', 'filled');




