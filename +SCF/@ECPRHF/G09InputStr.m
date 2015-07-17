function fileStr = G09InputStr(info)
newLine = sprintf('\n');
printCommand = sprintf('%s\n', ['#p uhf/genecp guess=harris', ...
    ' symmetry=none scf(NoVarAcc) iop(4/33=2) iop(3/33=1) iop(5/13=1)', ...
    ' scf(maxcycle=1)']);
printTitle = sprintf('%s\n', 'iop(3/33=1): print 1-e integrals; iop(5/13=1): do not terminate when scf fails');
printMol = sprintf('%d %d\n', info.chargeMult);
for iAtom = 1:size(info.cartesian, 1)
    printMol = [printMol, sprintf('%d %35.25f %35.25f %35.25f\n', info.cartesian(iAtom, :))]; %#ok
end
ecp = [fileread(info.ecpFile), newLine];
fileStr = [printCommand, newLine, printTitle, newLine, printMol, newLine, ecp, newLine, newLine];
end
