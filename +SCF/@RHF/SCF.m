function output = SCF(self, guess, diisType)
% initialize diis
cdiis20 = self.CDIIS(20);
ediis20 = self.EDIIS(20);
lciis20 = self.LCIIS(20);
cdiis6 = self.CDIIS(6);
ediis6 = self.EDIIS(6);
lciis6 = self.LCIIS(6);
adiis20 = self.ADIIS(20);

% initialize some variables
energy = 0;
maxErrSet = [];
energySet = [];
densVecSet = {};
for iter = 1:self.maxSCFIter
    tic
    oldEnergy = energy;
    
    if iter == 1
        if isfield(guess, 'guessOrbital')
            orbital = guess.guessOrbital;
            densVec = self.OrbToDensVec(orbital);
            fockVec = self.OrbToFockVec(orbital);
        elseif isfield(guess, 'guessDensVec')
            densVec = guess.guessDensVec;
            fockVec = self.DensVecToFockVec(densVec);
        else
            throw(MException('RHF:SCF', 'need an initial guess'));
        end
    else
        fockVec = self.OrbToFockVec(orbital);
    end
    energy = self.SCFEnergy(fockVec, densVec);
    
    % diis extrapolate Fock matrix
    cdiis20.Push(fockVec, densVec);
    ediis20.Push(fockVec, densVec, energy);
    lciis20.Push(fockVec, densVec);
    cdiis6.Push(fockVec, densVec);
    ediis6.Push(fockVec, densVec, energy);
    lciis6.Push(fockVec, densVec);
    adiis20.Push(fockVec, densVec);
    
    if iter == 2
        [~, ~, fockVecSet] = cdiis20.OptFockVector();
        fockVec = cdiis20.CalcFockVec([0.875; 0.125], fockVecSet);
    else
%     if true
        switch(diisType)
            case('ECmix20')
                [fockVec, maxErrSet] = ECmix(ediis20, cdiis20, cdiis20.MaxError(), maxErrSet);
            case('ECmix6')
                [fockVec, maxErrSet] = ECmix(ediis6, cdiis6, cdiis20.MaxError(), maxErrSet);
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
            case('LCe20')
                fockVec = ECe(lciis20, cdiis20, abs(energy - oldEnergy));
            case('L20')
                fockVec = lciis20.OptFockVector();
                disp('lciis(20)');
            case('L6')
                fockVec = lciis6.OptFockVector();
                disp('lciis(6)');
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
    end
    
    energySet(iter) = energy; %#ok
    fprintf('%0.8f\n',energy);
    
    oldDensVec = densVec;
    orbital = self.SolveFockVec(fockVec);
    densVec = self.OrbToDensVec(orbital);
    
    densVecSet{iter} = densVec; %#ok
    
    disp(mean(sqrt(mean((densVec - oldDensVec).^2))));
    disp(mean(max(abs(densVec - oldDensVec))));
    disp(abs(energy - oldEnergy));
    
    if mean(sqrt(mean((densVec - oldDensVec).^2))) < self.RMSDensityThreshold ...
            && mean(max(abs(densVec - oldDensVec))) < self.MaxDensityThreshold ...
            && abs(energy - oldEnergy) < self.EnergyThreshold
        break;
    end
    
    disp(['done iter ', num2str(iter)])
    toc
end

output.energy = energy;
output.energySet = energySet;
output.densVecSet = densVecSet;
output.iter = iter;
output.orbital = orbital;

end


function [fockVec, maxErrSet] = ECmix(ediis, cdiis, maxErr, maxErrSet)
maxErrSet = [maxErrSet, maxErr];
maxErr = maxErr * 2;
if maxErr > 1e-1 || maxErr > 1.1 * 2 * min(maxErrSet)
    fockVec = ediis.OptFockVector();
    disp(class(ediis));
elseif maxErr < 1e-4
    fockVec = cdiis.OptFockVector();
    disp(class(cdiis));
else
    [~, cdiisCoeffs, fockVecSet] = cdiis.OptFockVector();
    [~, ediisCoeffs, ~] = ediis.OptFockVector(length(cdiisCoeffs));
    coeffs = 10 .* maxErr .* ediisCoeffs + (1 - 10 .* maxErr) .* cdiisCoeffs;
    fockVec = cdiis.CalcFockVec(coeffs, fockVecSet);
    disp(['mix ', class(ediis), ' ', class(cdiis)]);
end
end

function [fockVec, maxErrSet] = EC(ediis, cdiis, maxErr, maxErrSet, iter)
if maxErr ~= 0
    maxErrSet(iter) = maxErr;
else
    maxErrSet(iter) = 1;
end
% if(maxErr > 1e-5 || maxErr > 1.1 * min(maxErrSet))
if maxErr > 1e-2
    fockVec = ediis.OptFockVector();
    disp(class(ediis));
else
    fockVec = cdiis.OptFockVector();
    disp(class(cdiis));
end
end

function fockVec = ECe(ediis, cdiis, energyDiff)
if energyDiff > 1e-3
    fockVec = ediis.OptFockVector();
    disp(class(ediis));
else
    fockVec = cdiis.OptFockVector();
    disp(class(cdiis));
end
end


