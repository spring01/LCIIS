classdef UKS < SCF.UHF & SCF.RKS
    
    methods
        
        function self = UKS(info)
            self@SCF.UHF(info);
            self@SCF.RKS(info);
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(self, orbital)
            occOrbA = orbital{1}(:, 1:self.numElectrons(1));
            occOrbB = orbital{2}(:, 1:self.numElectrons(2));
            self.currentV = self.matpsi.DFT_OccOrbToV(occOrbA, occOrbB);
            
            if self.hfExcCoeff == 0 % pure DFT, do not need K
                jMat = self.matpsi.JK_OccOrbToJ(occOrbA, occOrbB);
                gMat = repmat(sum(jMat, 3), [1 1 2]) + self.currentV;
            else
                self.matpsi.JK_CalcAllFromOccOrb(occOrbA, occOrbB);
                jMat = self.matpsi.JK_RetrieveJ();
                kMat = self.matpsi.JK_RetrieveK();
                gMat = repmat(sum(jMat, 3), [1 1 2]) + self.currentV ...
                     - self.hfExcCoeff * kMat;
            end
            
            fockVec = repmat(reshape(self.coreHamilt, [], 1), 1, 2) + ...
                [reshape(gMat(:,:,1), [], 1), reshape(gMat(:,:,2), [], 1)];
        end
        
        function fockVec = DensVecToFockVec(self, densVec)
            densMatA = reshape(densVec(:, 1), size(self.overlapMat));
            densMatB = reshape(densVec(:, 2), size(self.overlapMat));
            self.currentV = self.matpsi.DFT_DensToV(densMatA, densMatB);
            
            if self.hfExcCoeff == 0 % pure DFT, do not need K
                jMat = self.matpsi.JK_DensToJ(densMatA, densMatB);
                gMat = repmat(sum(jMat, 3), [1 1 2]) + self.currentV;
            else
                self.matpsi.JK_CalcAllFromOccOrb(densMatA, densMatB);
                jMat = self.matpsi.JK_RetrieveJ();
                gMat = repmat(sum(jMat, 3), [1, 1, 2]) + self.currentV ...
                     - self.hfExcCoeff .* self.matpsi.JK_RetrieveK();
            end
            
            fockVec = repmat(reshape(self.coreHamilt, [], 1), 1, 2) + ...
                [reshape(gMat(:,:,1), [], 1), reshape(gMat(:,:,2), [], 1)];
        end
        
        function energy = SCFEnergy(self, fockVec, densVec)
            energy = self.SCFEnergy@SCF.UHF(fockVec, densVec) ...
                - reshape(self.currentV(:,:,1), 1, []) * densVec(:, 1) / 2 ...
                - reshape(self.currentV(:,:,2), 1, []) * densVec(:, 2) / 2 ...
                + self.matpsi.DFT_EnergyXC();
        end
        
    end
    
end