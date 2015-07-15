classdef RKS < SCF.RHF
    
    properties (Access = protected)
        
        currentV;
        hfExcCoeff;
        
    end
    
    methods
        
        function obj = RKS(properties, dft)
            obj = obj@SCF.RHF(properties);
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
            obj.currentV = obj.matpsi2.DFT_OccOrbToV(occOrb);            
            if(obj.hfExcCoeff == 0) % pure DFT, do not need K
                gMat = 2 .* obj.matpsi2.JK_OccOrbToJ(occOrb) + obj.currentV;
            else % hybrid, do need K
                obj.matpsi2.JK_CalcAllFromOccOrb(occOrb);
                gMat = 2 .* obj.matpsi2.JK_RetrieveJ() + obj.currentV ...
                    - obj.hfExcCoeff * obj.matpsi2.JK_RetrieveK();
            end
            fockVec = reshape(obj.coreHamilt, [], 1) + reshape(gMat, [], 1);
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = obj.SCFEnergy@SCF.RHF(fockVec, densVec) ...
                - reshape(obj.currentV, 1, []) * densVec ...
                + obj.matpsi2.DFT_EnergyXC();
        end
        
    end
    
end