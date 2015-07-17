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
            inv_S_Half = inv(sqrtm(obj.overlapMat));
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
            obj.matpsi2.JK_CalcAllFromOccOrb(occOrb);
            gMat = 2 .* obj.matpsi2.JK_RetrieveJ() - obj.matpsi2.JK_RetrieveK();
            fockVec = reshape(obj.coreHamilt, [], 1) + reshape(gMat, [], 1);
        end
        
        function [orbital, orbEigValues] = SolveFockVec(~, fockVec, inv_S_Half)
            fockMat = reshape(fockVec, sqrt(numel(fockVec)), []);
            fockMat = inv_S_Half' * fockMat * inv_S_Half;
            fockMat = (fockMat + fockMat') ./ 2;
            [orbitalOtho, orbEigValues] = eig(fockMat);
            [orbEigValues, ascend_order] = sort(diag(orbEigValues));
            orbitalOtho = orbitalOtho(:, ascend_order);
            orbital = inv_S_Half * orbitalOtho;
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = (reshape(obj.coreHamilt, [], 1) + fockVec)'*densVec + obj.nucRepEnergy;
        end
        
        function cdiis = CDIIS(obj, numVectors)
            cdiis = SCF.CDIIS(obj.overlapMat, numVectors, 'r');
        end
        
        function ediis = EDIIS(obj, numVectors)
            ediis = SCF.EDIIS(obj.coreHamilt, numVectors, 'r');
        end
        
        function adiis = ADIIS(obj, numVectors)
            adiis = SCF.ADIIS(obj.coreHamilt, numVectors, 'r');
        end
        
        function lciis = LCIIS(obj, numVectors)
            lciis = SCF.LCIIS(obj.overlapMat, numVectors, 'r');
        end
        
    end
    
    methods (Static)
                
        function properties = MatPsi2Interface(matpsi2)
            properties.overlapMat = matpsi2.Integrals_Overlap();
            properties.coreHamilt = matpsi2.Integrals_Kinetic() + matpsi2.Integrals_Potential();
            properties.nucRepEnergy = matpsi2.Molecule_NucRepEnergy();
            properties.matpsi2 = matpsi2;
            chargeMult = matpsi2.Molecule_ChargeMult();
            properties.numElectrons = SCF.RHF.CalcNumElectrons( ...
                matpsi2.Molecule_NumElectrons(), chargeMult(2));
        end
        
    end
    
    methods(Static, Access = protected)
        
        function numElectrons = CalcNumElectrons(numTotalElectrons, mult)
            numAlphaElectrons = (numTotalElectrons + mult - 1) / 2;
            if(rem(numAlphaElectrons, 1) ~= 0)
                throw(MException('RHF:MatPsi2Interface', 'Number of electrons and multiplicity do not agree'));
            end
            numBetaElectrons = numTotalElectrons - numAlphaElectrons;
            numElectrons = [numAlphaElectrons; numBetaElectrons];
        end
        
    end
    
end