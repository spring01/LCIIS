function [energy, energySet, iter, orbital] = SCF(obj, guessOrbital, diisType)
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);
orbital = guessOrbital;
densVec = obj.OrbToDensVec(orbital);

% diis
cdiis20 = obj.CDIIS(20);
ediis20 = obj.EDIIS(20);
lciis20 = obj.LCIIS(20);
cdiis6 = obj.CDIIS(6);
ediis6 = obj.EDIIS(6);
lciis6 = obj.LCIIS(6);
adiis20 = obj.ADIIS(20);

% initialize some variables
energy = 0;
maxErrSet = [];
energySet = [];
for iter = 1:obj.maxSCFIter
    oldEnergy = energy;
    fockVec = obj.OrbToFockVec(orbital);
    energy = obj.SCFEnergy(fockVec, densVec);
    
    % diis extrapolate Fock matrix
    cdiis20.Push(fockVec, densVec);
    ediis20.Push(fockVec, densVec, energy);
    lciis20.Push(fockVec, densVec);
    cdiis6.Push(fockVec, densVec);
    ediis6.Push(fockVec, densVec, energy);
    lciis6.Push(fockVec, densVec);
    adiis20.Push(fockVec, densVec);
    switch(diisType)
        case('ECmix20')
            [fockVec, maxErrSet] = ECmix(ediis20, cdiis20, cdiis20.MaxError(), maxErrSet, iter);
        case('ECmix6')
            [fockVec, maxErrSet] = ECmix(ediis6, cdiis6, cdiis20.MaxError(), maxErrSet, iter);
        case('C20')
            fockVec = cdiis20.OptFockVector();
            disp('cdiis(20)');
        case('C6')
            fockVec = cdiis6.OptFockVector();
            disp('cdiis(6)');
        case('E20')
            fockVec = ediis20.OptFockVector();
            disp('ediis(20)')
        case('E6')
            fockVec = ediis6.OptFockVector();
            disp('ediis(6)')
        case('EC20')
            [fockVec, maxErrSet] = EC(ediis20, cdiis20, cdiis20.MaxError(), maxErrSet, iter);
        case('EC6')
            [fockVec, maxErrSet] = EC(ediis6, cdiis6, cdiis20.MaxError(), maxErrSet, iter);
        case('ECe6')
            fockVec = ECe(ediis6, cdiis6, abs(energy - oldEnergy));
        case('ECe20')
            fockVec = ECe(ediis20, cdiis20, abs(energy - oldEnergy));
        case('L20')
            fockVec = lciis20.OptFockVector();
            disp('lciis(20)');
        case('L6')
            fockVec = lciis6.OptFockVector();
            disp('mciis(6)');
        case('EL20')
            [fockVec, maxErrSet] = EC(ediis20, lciis20, cdiis20.MaxError(), maxErrSet, iter);
        case('EL6')
            [fockVec, maxErrSet] = EC(ediis6, lciis6, cdiis20.MaxError(), maxErrSet, iter);
        case('ELe6')
            fockVec = ECe(ediis6, lciis6, abs(energy - oldEnergy));
        case('ELe20')
            fockVec = ECe(ediis20, lciis20, abs(energy - oldEnergy));
        case('A20')
            fockVec = adiis20.OptFockVector();
            disp('adiis(20)')
        case('AC20')
            [fockVec, maxErrSet] = EC(adiis20, cdiis20, cdiis20.MaxError(), maxErrSet, iter);
        case('ACe20')
            fockVec = ECe(adiis20, cdiis20, abs(energy - oldEnergy));
        case('AL20')
            [fockVec, maxErrSet] = EC(adiis20, lciis20, cdiis20.MaxError(), maxErrSet, iter);
        case('ALe20')
            fockVec = ECe(adiis20, lciis20, abs(energy - oldEnergy));
    end
    
    energySet(iter) = energy; %#ok
    fprintf('%0.8f\n',energy);
    
    oldDensVec = densVec;
    orbital = obj.SolveFockVec(fockVec, inv_S_Half);
    densVec = obj.OrbToDensVec(orbital);
    
    disp(mean(sqrt(mean((densVec - oldDensVec).^2))));
    disp(mean(max(abs(densVec - oldDensVec))));
    disp(abs(energy - oldEnergy));
    
    if(mean(sqrt(mean((densVec - oldDensVec).^2))) < obj.RMSDensityThreshold ...
            && mean(max(abs(densVec - oldDensVec))) < obj.MaxDensityThreshold ...
            && abs(energy - oldEnergy) < obj.EnergyThreshold)
        break;
    end
    
    disp(['done iter ', num2str(iter)])
end

end


function [fockVec, maxErrSet] = ECmix(ediis, cdiis, maxErr, maxErrSet, iter)
if(maxErr ~= 0)
    maxErrSet(iter) = maxErr;
else
    maxErrSet(iter) = 1;
end
if(maxErr > 1e-1 || maxErr > 1.1 * min(maxErrSet))
    fockVec = ediis.OptFockVector();
    disp(class(ediis));
elseif(maxErr < 1e-4)
    fockVec = cdiis.OptFockVector();
    disp(class(cdiis));
else
    [~, ediisCoeffs, ~] = ediis.OptFockVector();
    [~, cdiisCoeffs, fockVecSet] = cdiis.OptFockVector();
    coeffs = 10.*maxErr .* ediisCoeffs + (1 - 10.*maxErr) .* cdiisCoeffs;
    fockVec = cdiis.CalcFockVec(coeffs, fockVecSet);
    disp(['mix ', class(ediis), ' ', class(cdiis)]);
end
end

function [fockVec, maxErrSet] = EC(ediis, cdiis, maxErr, maxErrSet, iter)
if(maxErr ~= 0)
    maxErrSet(iter) = maxErr;
else
    maxErrSet(iter) = 1;
end
% if(maxErr > 1e-5 || maxErr > 1.1 * min(maxErrSet))
if(maxErr > 1e-2)
    fockVec = ediis.OptFockVector();
    disp(class(ediis));
else
    fockVec = cdiis.OptFockVector();
    disp(class(cdiis));
end
end

function fockVec = ECe(ediis, cdiis, energyDiff)
if(energyDiff > 1e-3)
    fockVec = ediis.OptFockVector();
    disp(class(ediis));
else
    fockVec = cdiis.OptFockVector();
    disp(class(cdiis));
end
end


