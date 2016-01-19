function fileStr = G09InputStr(self)
% Only used to generate Harris guess and ECP matrix
% info.charge
% info.multiplicity
% info.cartesian
% info.g09BasisSetKeyword (optinal; genecp if ecp)
% info.ecpFile (optional)
info = self.info;
if ~isfield(info, 'g09BasisSetKeyword')
    info.g09BasisSetKeyword = info.basisSet;
end
newLine = sprintf('\n');
command = sprintf('%s\n', ...
    ['#p uhf/', info.g09BasisSetKeyword, ...
    ' guess=harris symmetry=none scf(NoVarAcc) iop(4/33=2) iop(3/33=1)', ...
    ' iop(5/13=1) scf(maxcycle=1)']);
title = sprintf('%s\n', ...
    ['iop(3/33=1): print 1-e integrals; ', ...
    'iop(5/13=1): do not terminate when scf fails']);
molecule = sprintf('%d %d\n', [info.charge, info.multiplicity]);
for iAtom = 1:size(info.cartesian, 1)
    molecule = [molecule, ...
        sprintf('%d %35.25f %35.25f %35.25f\n', info.cartesian(iAtom, :))]; %#ok
end
if ~isfield(info, 'ecpFile')
    ecp = [];
else
    ecp = [fileread(info.ecpFile), newLine];
end
fileStr = [command, newLine, title, newLine, molecule, newLine, ...
    ecp, newLine, newLine];
end
