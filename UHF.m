classdef UHF < RHF
    
    methods
        
        function obj = UHF(properties)
            obj = obj@RHF(properties);
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
            jMat = obj.matpsi2.JK_OccOrbToJ(occOrbAlpha, occOrbBeta);
            kMat = obj.matpsi2.JK_OccOrbToK(occOrbAlpha, occOrbBeta);
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
                    = obj.SolveFockVec@RHF(fockVec(:, spin), inv_S_Half);
            end
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = (reshape(obj.coreHamilt, [], 1)' * sum(densVec, 2) ...
                + fockVec(:, 1)' * densVec(:, 1) ...
                + fockVec(:, 2)' * densVec(:, 2)) / 2 + obj.nucRepEnergy;
        end
        
        function newVec = Damping(~, dampingCoeff, vec, oldVec)
            newVec = zeros(size(vec));
            coeffs = [dampingCoeff; (1 - dampingCoeff)];
            for spin = 1:2
                newVec(:, spin) = [vec(:, spin), oldVec(:, spin)] * coeffs;
            end
        end
        
        function cdiis = CDIIS(obj, numVectors)
            cdiis = CDIIS(obj.overlapMat, numVectors, 'u');
        end
        
        function ediis = EDIIS(obj, numVectors)
            ediis = EDIIS(obj.coreHamilt, numVectors, 'u');
        end
        
        function adiis = ADIIS(obj, numVectors)
            adiis = ADIIS(obj.coreHamilt, numVectors, 'u');
        end
        
        function lciis = LCIIS(obj, numVectors)
            lciis = LCIIS(obj.overlapMat, numVectors, 'u');
        end
        
    end
    
end