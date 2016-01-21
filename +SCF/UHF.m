classdef UHF < SCF.RHF
    
    methods
        
        function self = UHF(info)
            self = self@SCF.RHF(info);
        end
        
        function [densVec, orbital] = CoreGuess(self)
            coreFockVec = reshape(self.coreHamilt, [], 1);
            orbital = self.SolveFockVec([coreFockVec, coreFockVec]);
            densVec = self.OrbToDensVec(orbital);
        end
        
        function [densVec, orbital] = HarrisGuess(self)
            orbital = {self.G09ReadMatrix('HarrisGuessMOAlpha'), ...
                       self.G09ReadMatrix('HarrisGuessMOBeta')};
            [g09Order, g09Renormby] = self.G09OrderG09RenormBy();
            g09RenormByMat = repmat(g09Renormby, [1, size(orbital{1}, 2)]);
            orbital{1} = orbital{1} .* g09RenormByMat;
            orbital{1} = orbital{1}(g09Order, :);
            orbital{2} = orbital{2} .* g09RenormByMat;
            orbital{2} = orbital{2}(g09Order, :);
            densVec = self.OrbToDensVec(orbital);
        end
        
    end
    
    methods (Access = protected)
        
        function densVec = OrbToDensVec(self, orbital)
            densVec = zeros(numel(self.overlapMat), 2);
            occOrbA = orbital{1}(:, 1:self.numElectrons(1));
            densVec(:, 1) = reshape(occOrbA * occOrbA', [], 1);
            occOrbB = orbital{2}(:, 1:self.numElectrons(2));
            densVec(:, 2) = reshape(occOrbB * occOrbB', [], 1);
        end
        
        function fockVec = OrbToFockVec(self, orbital)
            occOrbAlpha = orbital{1}(:, 1:self.numElectrons(1));
            occOrbBeta = orbital{2}(:, 1:self.numElectrons(2));
            self.matpsi.JK_CalcAllFromOccOrb(occOrbAlpha, occOrbBeta);
            fockVec = self.RetrieveFockVec();
        end
        
        function fockVec = DensVecToFockVec(self, densVec)
            densMatA = reshape(densVec(:, 1), size(self.overlapMat));
            densMatB = reshape(densVec(:, 2), size(self.overlapMat));
            self.matpsi.JK_CalcAllFromDens(densMatA, densMatB);
            fockVec = self.RetrieveFockVec();
        end
        
        function [orbital, orbEigValues] = SolveFockVec(self, fockVec)
            [orbital{1}, orbEigValues{1}] ...
                = self.SolveFockVec@SCF.RHF(fockVec(:, 1));
            [orbital{2}, orbEigValues{2}] ...
                = self.SolveFockVec@SCF.RHF(fockVec(:, 2));
        end
        
        function energy = SCFEnergy(self, fockVec, densVec)
            energy = (reshape(self.coreHamilt, [], 1)' * sum(densVec, 2) ...
                + fockVec(:, 1)' * densVec(:, 1) ...
                + fockVec(:, 2)' * densVec(:, 2)) / 2 + self.nucRepEnergy;
        end
        
        function cdiis = CDIIS(self, numVectors)
            cdiis = SCF.CDIIS(self.overlapMat, numVectors, 'u');
        end
        
        function ediis = EDIIS(self, numVectors)
            ediis = SCF.EDIIS(self.coreHamilt, numVectors, 'u');
        end
        
        function adiis = ADIIS(self, numVectors)
            adiis = SCF.ADIIS(self.coreHamilt, numVectors, 'u');
        end
        
        function lciis = LCIIS(self, numVectors)
            lciis = SCF.ULCIIS(self.overlapMat, numVectors);
        end
        
    end
    
    methods (Access = private)
        
        function fockVec = RetrieveFockVec(self)
            jMat = self.matpsi.JK_RetrieveJ();
            gMat = repmat(sum(jMat, 3), [1, 1, 2]) - self.matpsi.JK_RetrieveK();
            fockVec = repmat(reshape(self.coreHamilt, [], 1), 1, 2) + ...
                [reshape(gMat(:, :, 1), [], 1), reshape(gMat(:, :, 2), [], 1)];
        end
        
    end
    
end