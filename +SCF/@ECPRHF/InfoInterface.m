function properties = InfoInterface(info)
matpsi2 = info.matpsi2;

info.chargeMult = matpsi2.Molecule_ChargeMult();
SCF.ECPRHF.RunG09(info);

order = SCF.ECPRHF.G09ToPsi4BasisOrder(matpsi2.BasisSet_ShellNumFunctions());
ecpMat = SCF.ECPRHF.G09ReadMatrix('ecpInt');
ecpMat = ecpMat(order, order);

properties.overlapMat = matpsi2.Integrals_Overlap();
screenedZxyz = [ScreenedAtomNumbers(matpsi2.Molecule_AtomicNumbers())', matpsi2.Molecule_Geometry()];
properties.coreHamilt = matpsi2.Integrals_Kinetic() ...
    + matpsi2.Integrals_PotentialPtQ(screenedZxyz) ...
    + ecpMat;
chargeMult = matpsi2.Molecule_ChargeMult();
properties.numElectrons = SCF.RHF.CalcNumElectrons( ...
    sum(screenedZxyz(:, 1)) - chargeMult(1), chargeMult(2));
properties.nucRepEnergy = NucRepEnergy(screenedZxyz);
properties.matpsi2 = matpsi2;
end


function atomNumbers = ScreenedAtomNumbers(atomNumbers)
for i = 1:length(atomNumbers)
    numCoreElectrons = 0;
    if(11 <= atomNumbers(i) && atomNumbers(i) <= 29)
        numCoreElectrons = 10;
    elseif(atomNumbers(i) == 30)
        numCoreElectrons = 18;
    elseif(31 <= atomNumbers(i) && atomNumbers(i) <= 47)
        numCoreElectrons = 28;
    elseif(atomNumbers(i) == 48)
        numCoreElectrons = 36;
    elseif(49 <= atomNumbers(i) && atomNumbers(i) <= 56)
        numCoreElectrons = 46;
    elseif(atomNumbers(i) == 92)
        numCoreElectrons = 78;
    end
    atomNumbers(i) = atomNumbers(i) - numCoreElectrons;
end
end

function nucRepEnergy = NucRepEnergy(cartesian)
distMat = dist(cartesian(:, 2:end)');
chargeMat = cartesian(:, 1) * cartesian(:, 1)';
nucRepEnergy = sum(chargeMat(tril(true(size(chargeMat)), -1)) ...
    ./ distMat(tril(true(size(distMat)), -1)));
end

