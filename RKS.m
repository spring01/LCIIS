classdef RKS < RHF
    
    properties (Access = protected)
        
        previousV;
        currentV;
        hfExcCoeff;
        
    end
    
    methods
        
        function obj = RKS(properties, dft)
            obj = obj@RHF(properties);
            obj.matpsi2.DFT_Initialize(dft);
            if(strcmpi(dft, 'b3lyp') || strcmpi(dft, 'b3lypv5'))
                obj.hfExcCoeff = 0.2;
            else
                obj.hfExcCoeff = 0;
            end
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            occOrb = orbital(:, 1:obj.numElectrons(1));
            obj.previousV = obj.currentV;
            obj.currentV = obj.matpsi2.DFT_OccOrbToV(occOrb);
            gMat = 2 .* obj.matpsi2.JK_OccOrbToJ(occOrb) + obj.currentV;
            if(obj.hfExcCoeff ~= 0)
                gMat = gMat - obj.hfExcCoeff * obj.matpsi2.JK_OccOrbToK(occOrb);
            end
            fockVec = reshape(obj.coreHamilt, [], 1) + reshape(gMat, [], 1);
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = obj.SCFEnergy@RHF(fockVec, densVec) ...
                - reshape(obj.currentV, 1, []) * densVec ...
                + obj.matpsi2.DFT_EnergyXC();
        end
        
        function energy = DampedSCFEnergy(obj, fockVec, densVec, dampingCoeff, ~)
            obj.matpsi2.DFT_DensToV(reshape(densVec, sqrt(numel(densVec)), []));
            dampedV = obj.currentV .* dampingCoeff + obj.previousV .* (1 - dampingCoeff);
            energy = (reshape(obj.coreHamilt, [], 1) + fockVec)'*densVec + obj.nucRepEnergy ...
                - reshape(dampedV, 1, []) * densVec + obj.matpsi2.DFT_EnergyXC();
        end
        
    end
    
end