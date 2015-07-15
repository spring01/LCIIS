classdef UHF < SCF.RHF
    
    methods
        
        function obj = UHF(properties)
            obj = obj@SCF.RHF(properties);
        end
        
    end
    
    methods (Access = protected)
        
        function densVec = OrbToDensVec(obj, orbital)
            densVec = zeros(numel(obj.overlapMat), 2);
            for spin = 1:2
                occOrb = orbital{spin}(:, 1:obj.numElectrons(spin));
                densVec(:, spin) = reshape(occOrb * occOrb', [], 1);
            end
        end
        
        function fockVec = OrbToFockVec(obj, orbital)
            occOrbAlpha = orbital{1}(:, 1:obj.numElectrons(1));
            occOrbBeta = orbital{2}(:, 1:obj.numElectrons(2));
            obj.matpsi2.JK_CalcAllFromOrb(occOrbAlpha, occOrbBeta);
            jMat = obj.matpsi2.JK_RetrieveJ();
            kMat = obj.matpsi2.JK_RetrieveK();
            gMat = repmat(sum(jMat, 3), [1 1 2]) - kMat;
            oeiVec = reshape(obj.coreHamilt, [], 1);
            fockVec = repmat(oeiVec, 1, 2) + ...
                [reshape(gMat(:,:,1), [], 1), reshape(gMat(:,:,2), [], 1)];
        end
        
        function [orbital, orbEigValues] = SolveFockVec(obj, fockVec, inv_S_Half)
            % when fockVec is a single vector, repeat it; otherwise "has no effect"
            fockVec = repmat(fockVec, 1, 2);
            
            orbEigValues = cell(1, 2);
            orbital = cell(1, 2);
            for spin = 1:2
                [orbital{spin}, orbEigValues{spin}] ...
                    = obj.SolveFockVec@SCF.RHF(fockVec(:, spin), inv_S_Half);
            end
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = (reshape(obj.coreHamilt, [], 1)' * sum(densVec, 2) ...
                + fockVec(:, 1)' * densVec(:, 1) ...
                + fockVec(:, 2)' * densVec(:, 2)) / 2 + obj.nucRepEnergy;
        end
        
        function cdiis = CDIIS(obj, numVectors)
            cdiis = SCF.CDIIS(obj.overlapMat, numVectors, 'u');
        end
        
        function ediis = EDIIS(obj, numVectors)
            ediis = SCF.EDIIS(obj.coreHamilt, numVectors, 'u');
        end
        
        function adiis = ADIIS(obj, numVectors)
            adiis = SCF.ADIIS(obj.coreHamilt, numVectors, 'u');
        end
        
        function lciis = LCIIS(obj, numVectors)
            lciis = SCF.LCIIS(obj.overlapMat, numVectors, 'u');
        end
        
    end
    
end