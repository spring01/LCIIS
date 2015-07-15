classdef RHF < handle
    
    properties (Access = protected)
        
        overlapMat;
        coreHamilt;
        nucRepEnergy;
        numElectrons;
        matpsi2;
        
        maxSCFIter = 200;
        RMSDensityThreshold = 1e-8;
        MaxDensityThreshold = 1e-6;
        EnergyThreshold = 1e-6;
        
    end
    
    methods
        
        function obj = RHF(properties)
            obj.overlapMat = properties.overlapMat;
            obj.coreHamilt = properties.coreHamilt;
            obj.nucRepEnergy = properties.nucRepEnergy;
            obj.numElectrons = properties.numElectrons;
            obj.matpsi2 = properties.matpsi2;
        end
        
        function [densVec, orbital] = CoreGuess(obj)
            coreFockVec = reshape(obj.coreHamilt, [], 1);
            inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);
            orbital = obj.SolveFockVec(coreFockVec, inv_S_Half);
            densVec = obj.OrbToDensVec(orbital);
        end
        
    end
    
    methods (Access = protected)
        
        function densVec = OrbToDensVec(obj, orbital)
            occOrb = orbital(:, 1:obj.numElectrons(1));
            densVec = reshape(occOrb * occOrb', [], 1);
        end
        
        function fockVec = OrbToFockVec(obj, orbital)
            occOrb = orbital(:, 1:obj.numElectrons(1));
            gMat = 2 .* obj.matpsi2.JK_OccOrbToJ(occOrb) ...
                - obj.matpsi2.JK_OccOrbToK(occOrb);
            fockVec = reshape(obj.coreHamilt, [], 1) + reshape(gMat, [], 1);
        end
        
        function [orbital, orbEigValues] = SolveFockVec(obj, fockVec, inv_S_Half)
            fockMat = reshape(fockVec, sqrt(numel(fockVec)), []);
            [orbitalOtho, orbEigValues] = eig(inv_S_Half*fockMat*inv_S_Half);
            [orbEigValues, ascend_order] = sort(diag(orbEigValues));
            orbitalOtho = obj.OrthDegenOrbitals(orbitalOtho(:, ascend_order), orbEigValues);
            orbital = inv_S_Half * orbitalOtho;
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = (reshape(obj.coreHamilt, [], 1) + fockVec)'*densVec + obj.nucRepEnergy;
        end
        
        function energy = DampedSCFEnergy(obj, fockVec, densVec, ~, ~)
            energy = obj.SCFEnergy(fockVec, densVec);
        end
        
        function newVec = Damping(~, dampingCoeff, vec, oldVec)
            coeffs = [dampingCoeff; (1 - dampingCoeff)];
            newVec = [vec, oldVec] * coeffs;
        end
        
        function cdiis = CDIIS(obj, numVectors)
            cdiis = CDIIS(obj.overlapMat, numVectors, 'r');
        end
        
        function ediis = EDIIS(obj, numVectors)
            ediis = EDIIS(obj.coreHamilt, numVectors, 'r');
        end
        
        function adiis = ADIIS(obj, numVectors)
            adiis = ADIIS(obj.coreHamilt, numVectors, 'r');
        end
        
        function lciis = LCIIS(obj, numVectors)
            lciis = LCIIS(obj.overlapMat, numVectors, 'r');
        end
        
    end
    
    methods (Access = private)
        
        function orbitalOtho = OrthDegenOrbitals(~, orbitalOtho, orbEigValues)
            diffVec = [0; diff(orbEigValues)];
            diffVec(abs(diffVec) < 1e-8) = 0;
            indVec = 1:length(orbEigValues);
            indStart = [1, indVec(diffVec~=0)];
            for iUniqEVal = 1:length(indStart) - 1
                indRange = indStart(iUniqEVal):indStart(iUniqEVal+1);
                temp = orbitalOtho(:, indRange);
                for iSubspace = 2:size(temp, 2)
                    for jSubspace = 1:iSubspace-1
                        temp(:, iSubspace) = temp(:, iSubspace) - (temp(:, iSubspace)'*temp(:, jSubspace)) .* temp(:, jSubspace);
                    end
                    temp(:, iSubspace) = temp(:, iSubspace) ./ norm(temp(:, iSubspace));
                end
                orbitalOtho(:, indRange) = temp;
            end
        end
        
    end
    
    methods (Static)
                
        function properties = MatPsi2Interface(matpsi2)
            properties.overlapMat = matpsi2.Integrals_Overlap();
            properties.coreHamilt = matpsi2.Integrals_Kinetic() + matpsi2.Integrals_Potential();
            properties.nucRepEnergy = matpsi2.Molecule_NucRepEnergy();
            properties.matpsi2 = matpsi2;
            
            chargeMult = matpsi2.Molecule_ChargeMult();
            mult = chargeMult(2);
            numTotalElectrons = matpsi2.Molecule_NumElectrons();
            numAlphaElectrons = (numTotalElectrons + mult - 1) / 2;
            if(rem(numAlphaElectrons, 1) ~= 0)
                throw(MException('RHF:MatPsi2Interface', 'Number of electrons and multiplicity do not agree'));
            end
            numBetaElectrons = numTotalElectrons - numAlphaElectrons;
            properties.numElectrons = [numAlphaElectrons; numBetaElectrons];
        end
        
    end
    
end