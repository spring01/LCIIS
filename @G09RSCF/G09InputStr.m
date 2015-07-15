function fileStr = G09InputStr(info)
printOrbAlpha = [];
guessStr = ' guess=core';
if(isfield(info, 'orbAlpha'))
    guessStr = ' guess=cards';
    printOrbAlpha = sprintf('%s\n', '(1d35.25)');
    for iOrb = 1:size(info.orbAlpha, 2)
        printOrbAlpha = [printOrbAlpha, sprintf('%d\n', iOrb)]; %#ok
        printOrbAlpha = [printOrbAlpha, sprintf('%35.25d\n', info.orbAlpha(:, iOrb))]; %#ok
    end
    printOrbAlpha = [printOrbAlpha, sprintf('%d\n', 0)];
elseif(isfield(info, 'harris'))
    guessStr = ' guess=harris iop(4/33=3)';
end
printOrbBeta = [];
if(isfield(info, 'orbBeta'))
    for iOrb = 1:size(info.orbBeta, 2)
        printOrbBeta = [printOrbBeta, sprintf('%d\n', iOrb)]; %#ok
        printOrbBeta = [printOrbBeta, sprintf('%35.25d\n', info.orbBeta(:, iOrb))]; %#ok
    end
end

maxCycle = '1';
dampingNum = 100;
if(isfield(info, 'dampingCoeff'))
    maxCycle = '2';
    dampingNum = info.dampingCoeff * 100;
end

newLine = sprintf('\n');
printCommand = sprintf('%s\n', ['#p ', info.method, '/', info.basisSet , guessStr, ...
    ' symmetry=none population=full scf(NoVarAcc) iop(5/33=3) iop(3/33=1) iop(5/13=1)', ...
    ' scf(maxcycle=', maxCycle, ') iop(5/18=', num2str(dampingNum), ') IOp(4/6=200000000) IOp(5/14=1)']);
printTitle = sprintf('%s\n', 'iop(5/33=3): print fock; iop(3/33=1): print 1-e integrals; iop(5/13=1): do not terminate when scf fails');
printMol = sprintf('%d %d\n', info.chargeMult);
for iAtom = 1:size(info.cartesian, 1)
    printMol = [printMol, sprintf('%d %35.25f %35.25f %35.25f\n', info.cartesian(iAtom, :))]; %#ok
end

ecp = [];
if(isfield(info, 'ecpFile'))
    ecp = [fileread(info.ecpFile), newLine];
end

fileStr = [printCommand, newLine, printTitle, newLine, printMol, newLine, ecp, printOrbAlpha, printOrbBeta, newLine, newLine];
end
