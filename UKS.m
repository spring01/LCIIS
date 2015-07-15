classdef UKS < UHF & RKS
    
    methods
        
        function obj = UKS(properties, dft)
            obj@RKS(properties, dft);
            obj@UHF(properties);
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            occOrbAlpha = orbital{1}(:, 1:obj.numElectrons(1));
            occOrbBeta = orbital{2}(:, 1:obj.numElectrons(2));
            jMat = obj.matpsi2.JK_OccOrbToJ(occOrbAlpha, occOrbBeta);
            obj.previousV = obj.currentV;
            obj.currentV = obj.matpsi2.DFT_OccOrbToV(occOrbAlpha, occOrbBeta);
            gMat = repmat(sum(jMat, 3), [1 1 2]) + obj.currentV;
            if(obj.hfExcCoeff ~= 0)
                kMat = obj.matpsi2.JK_OccOrbToK(occOrbAlpha, occOrbBeta);
                gMat = gMat - obj.hfExcCoeff * kMat;
            end
            oeiVec = reshape(obj.coreHamilt, [], 1);
            fockVec = repmat(oeiVec, 1, 2) + ...
                [reshape(gMat(:,:,1), [], 1), reshape(gMat(:,:,2), [], 1)];
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = obj.SCFEnergy@UHF(fockVec, densVec) ...
                - reshape(obj.currentV(:,:,1), 1, []) * densVec(:, 1) / 2 ...
                - reshape(obj.currentV(:,:,2), 1, []) * densVec(:, 2) / 2 ...
                + obj.matpsi2.DFT_EnergyXC();
        end
        
        function energy = DampedSCFEnergy(obj, fockVec, densVec, dampingCoeff, ~)
            densAlphaMat = reshape(densVec(:, 1), sqrt(numel(densVec(:, 1))), []);
            densBetaMat = reshape(densVec(:, 2), sqrt(numel(densVec(:, 2))), []);
            obj.matpsi2.DFT_DensToV(densAlphaMat, densBetaMat);
            dampedV = obj.currentV .* dampingCoeff + obj.previousV .* (1 - dampingCoeff);
            energy = (reshape(obj.coreHamilt, [], 1)' * sum(densVec, 2) ...
                + fockVec(:, 1)' * densVec(:, 1) ...
                + fockVec(:, 2)' * densVec(:, 2)) / 2 + obj.nucRepEnergy ...
                - reshape(dampedV(:,:,1), 1, []) * densVec(:, 1) / 2 ...
                - reshape(dampedV(:,:,2), 1, []) * densVec(:, 2) / 2 ...
                + obj.matpsi2.DFT_EnergyXC();
        end
        
    end
    
end