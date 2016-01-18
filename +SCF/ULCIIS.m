classdef ULCIIS < SCF.RLCIIS
    
    methods
        
        function self = ULCIIS(overlapMatrix, maxNumPairs)
            self@SCF.RLCIIS(overlapMatrix, maxNumPairs)
            self.comms = [self.comms; self.comms];
        end
        
    end
    
    methods (Access = protected)
        
        function comm = CommBetween(self, ind1, ind2)
            commA = self.Commutator(self.pairs{ind1}.fockVector(:, 1), ...
                                    self.pairs{ind2}.densVector(:, 1));
            commB = self.Commutator(self.pairs{ind1}.fockVector(:, 2), ...
                                    self.pairs{ind2}.densVector(:, 2));
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
