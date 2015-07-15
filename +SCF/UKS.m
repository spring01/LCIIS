classdef UKS < SCF.UHF & SCF.RKS
    
    methods
        
        function obj = UKS(properties, dft)
            obj@SCF.RKS(properties, dft);
            obj@SCF.UHF(properties);
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            occOrbAlpha = orbital{1}(:, 1:obj.numElectrons(1));
            occOrbBeta = orbital{2}(:, 1:obj.numElectrons(2));
            obj.currentV = obj.matpsi2.DFT_OccOrbToV(occOrbAlpha, occOrbBeta);
            
            if(obj.hfExcCoeff == 0) % pure DFT, do not need K
                jMat = obj.matpsi2.JK_OccOrbToJ(occOrbAlpha, occOrbBeta);
                gMat = repmat(sum(jMat, 3), [1 1 2]) + obj.currentV;
            else
                obj.matpsi2.JK_CalcAllFromOrb(occOrbAlpha, occOrbBeta);
                jMat = obj.matpsi2.JK_RetrieveJ();
                kMat = obj.matpsi2.JK_RetrieveK();
                gMat = repmat(sum(jMat, 3), [1 1 2]) + obj.currentV ...
                    - obj.hfExcCoeff * kMat;
            end
            
            oeiVec = reshape(obj.coreHamilt, [], 1);
            fockVec = repmat(oeiVec, 1, 2) + ...
                [reshape(gMat(:,:,1), [], 1), reshape(gMat(:,:,2), [], 1)];
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = obj.SCFEnergy@SCF.UHF(fockVec, densVec) ...
                - reshape(obj.currentV(:,:,1), 1, []) * densVec(:, 1) / 2 ...
                - reshape(obj.currentV(:,:,2), 1, []) * densVec(:, 2) / 2 ...
                + obj.matpsi2.DFT_EnergyXC();
        end
        
    end
    
end