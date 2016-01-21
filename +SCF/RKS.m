classdef RKS < SCF.RHF
    
    properties (Access = protected)
        
        currentV;
        hfExcCoeff;
        
    end
    
    methods
        
        function self = RKS(info)
            self = self@SCF.RHF(info);
            self.matpsi.DFT_Initialize(info.dft);
            if(strcmpi(info.dft, 'b3lyp') || strcmpi(info.dft, 'b3lypv5'))
                self.hfExcCoeff = 0.2;
            else
                self.hfExcCoeff = 0;
            end
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(self, orbital)
            occOrb = orbital(:, 1:self.numElectrons(1));
            self.currentV = self.matpsi.DFT_OccOrbToV(occOrb);            
            if self.hfExcCoeff == 0 % pure DFT, do not need K
                gMat = 2 .* self.matpsi.JK_OccOrbToJ(occOrb) + self.currentV;
            else % hybrid, do need K
                self.matpsi.JK_CalcAllFromOccOrb(occOrb);
                gMat = 2 .* self.matpsi.JK_RetrieveJ() + self.currentV ...
                     - self.hfExcCoeff .* self.matpsi.JK_RetrieveK();
            end
            fockVec = reshape(self.coreHamilt + gMat, [], 1);
        end
        
        function fockVec = DensVecToFockVec(self, densVec)
            densMat = reshape(densVec, size(self.overlapMat));
            self.currentV = self.matpsi.DFT_DensToV(densMat);            
            if self.hfExcCoeff == 0 % pure DFT, do not need K
                gMat = 2 .* self.matpsi.JK_DensToJ(densMat) + self.currentV;
            else % hybrid, do need K
                self.matpsi.JK_CalcAllFromDens(densMat);
                gMat = 2 .* self.matpsi.JK_RetrieveJ() + self.currentV ...
                     - self.hfExcCoeff .* self.matpsi.JK_RetrieveK();
            end
            fockVec = reshape(self.coreHamilt + gMat, [], 1);
        end
        
        function energy = SCFEnergy(self, fockVec, densVec)
            energy = self.SCFEnergy@SCF.RHF(fockVec, densVec) ...
                - reshape(self.currentV, 1, []) * densVec ...
                + self.matpsi.DFT_EnergyXC();
        end
        
    end
    
end