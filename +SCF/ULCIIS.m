classdef ULCIIS < SCF.RLCIIS
    
    methods
        
        function self = ULCIIS(overlapMatrix, maxNumPairs)
            self@SCF.RLCIIS(overlapMatrix, maxNumPairs)
            self.comms = [self.comms; self.comms];
        end
        
    end
    
    methods (Access = protected)
        
        function pair = WrapFockDensInPair(self, newFockVector, newDensVector)
            % This function defines interface.
            % In unrestricted SCF, newFockVector stores the Fock matrices, and
            % newDensVector stores the density matrices, both as two nbf^2 by 1
            % vectors; newFockVector(:, 1) is the alpha Fock matrix, and
            % newFockVector(:, 2) is the beta Fock matrix; so as
            % newDensVector.
            nbf = size(self.sqrtmS, 1);
            pair.fockVector = newFockVector;
            fockOrthoMat = cell(1, size(newDensVector, 2));
            for i = 1:size(newFockVector, 2)
                fockMat = reshape(newFockVector(:, i), nbf, nbf);
                fockOrthoMat{i} = self.invSqrtmS * fockMat * self.invSqrtmS;
            end
            densOrthoMat = cell(1, size(newDensVector, 2));
            for i = 1:size(newDensVector, 2)
                densMat = reshape(newDensVector(:, i), nbf, nbf);
                densOrthoMat{i} = self.sqrtmS * densMat * self.sqrtmS;
            end
            pair.fockOrthoMatrices = fockOrthoMat;
            pair.densOrthoMatrices = densOrthoMat;
        end
        
        function comm = CommBetween(self, ind1, ind2)
            commA = self.Commutator(self.pairs{ind1}.fockOrthoMatrices{1}, ...
                                    self.pairs{ind2}.densOrthoMatrices{1});
            commB = self.Commutator(self.pairs{ind1}.fockOrthoMatrices{2}, ...
                                    self.pairs{ind2}.densOrthoMatrices{2});
            comm = [commA; commB];
        end
        
        function [newFockVector, fockVectors] = ExtrapolateFock(self, coeffVec)
            numPairs = length(coeffVec);
            start = length(self.pairs) - numPairs;
            lenFockVector = size(self.pairs{end}.fockVector, 1);
            fockVectors{1} = zeros(lenFockVector, numPairs);
            fockVectors{2} = zeros(lenFockVector, numPairs);
            for i = 1:numPairs
                fockVectors{1}(:, i) = self.pairs{start + i}.fockVector(:, 1);
                fockVectors{2}(:, i) = self.pairs{start + i}.fockVector(:, 2);
            end
            newFockVector = zeros(lenFockVector, 2);
            newFockVector(:, 1) = fockVectors{1} * coeffVec;
            newFockVector(:, 2) = fockVectors{2} * coeffVec;
        end
        
    end
end
