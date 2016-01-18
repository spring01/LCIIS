classdef RLCIIS < handle
    
    properties (Access = protected)
        
        pairs;
        commutators;
        
    end
    
    properties (Access = private)
        
        maxNumPairs;
        
        pairsFull;
        bigMatrix;
        tensor;
        overlap;
        sqrtmOverlap;
        invSqrtmOverlap;
        
        tempTensorForGrad;
        tempTensorForHess;
        
        maxIterNewton = 200;
        minHessRcond = 1e-25;
        gradNormThres = 1e-12;
        
    end
    
    methods
        
        function self = RLCIIS(overlapMatrix, maxNumPairs)
            if nargin < 2
                maxNumPairs = 5;
            end
            self.pairs = {};
            self.maxNumPairs = maxNumPairs;
            self.commutators = zeros(numel(overlapMatrix), maxNumPairs^2);
            self.bigMatrix = zeros(maxNumPairs^2, maxNumPairs^2);
            self.sqrtmOverlap = sqrtm(overlapMatrix);
            self.invSqrtmOverlap = inv(self.sqrtmOverlap);
        end
        
        function Push(self, newFockVector, newDensVector)
            newPair.fockVector = newFockVector;
            newPair.densVector = newDensVector;
            if length(self.pairs) == size(self.commutators, 1)
                self.pairs(1:end-1) = self.pairs(2:end);
                self.pairs{end} = newPair;
                self.pairsFull = true;
            else
                self.pairs{end+1} = newPair;
                self.pairsFull = false;
            end
        end
        
        function [optFockVector, coeffVec, useFockVectors] = OptFockVector(self)
            if self.pairsFull
                self.MoveCommutatorsFull();
            else
                self.MoveCommutatorsNotFull();
            end
            self.UpdateCommutators();
            self.CalcTensor();
            coeffVec = self.FindCoeffVec();
            [optFockVector, useFockVectors] = self.ExtrapolateFock(coeffVec);
        end
        
    end
    
    methods (Access = protected)
        
        function commutator = CommutatorForInd(self, ind1, ind2)
            commutator = self.Commutator(self.pairs{ind1}.fockVector, ...
                                         self.pairs{ind2}.densVector);
        end
        
        function comm = Commutator(self, singleFockVector, singleDensVector)
            nbf = size(self.sqrtmOverlap, 1);
            fockMat = reshape(singleFockVector, nbf, nbf);
            densMat = reshape(singleDensVector, nbf, nbf);
            fockDens = self.invSqrtmOverlap * fockMat ...
                     * densMat * self.sqrtmOverlap;
            comm = reshape(fockDens - fockDens', [], 1);
        end
        
        function [newFockVector, fockVectors] = ExtrapolateFock(self, coeffVec)
            numPairs = length(coeffVec);
            start = length(self.pairs) - numPairs;
            fockVectors = zeros(length(self.pairs{end}.fockVector), numPairs);
            for i = 1:numPairs
                fockVectors(:, i) = self.pairs{start + i}.fockVector;
            end
            newFockVector = fockVectors * coeffVec;
        end
        
    end
    
    methods (Access = private)
        
        function CalcTensor(self)
            numPairs = length(self.pairs);
            comms = self.commutators(:, 1:numPairs^2);
            self.tensor = reshape(comms' * comms, numPairs * ones(1, 4));
        end
        
        function MoveCommutatorsFull(self)
            numPairs = self.maxNumPairs;
            for i = 2:numPairs
                srcInd = ((i - 1) * numPairs + 2):(i * numPairs);
                tgtInd = srcInd - numPairs - 1;
                self.commutators(:, tgtInd) = self.commutators(:, srcInd);
            end
        end
        
        function MoveCommutatorsNotFull(self)
            numPairs = length(self.pairs);
            for i = (numPairs - 1):-1:2
                srcInd = ((i - 1) * (numPairs - 1) + 1):(i * (numPairs - 1));
                tgtInd = srcInd + i - 1;
                self.commutators(:, tgtInd) = self.commutators(:, srcInd);
            end
        end
        
        function UpdateCommutators(self)
            numPairs = length(self.pairs);
            for i = 1:(numPairs-1)
                self.commutators(:, i * numPairs) = ...
                    self.CommutatorForInd(i, numPairs);
            end
            
            for ind = ((numPairs - 1) * numPairs + 1):numPairs^2
                self.commutators(:, ind) = ...
                    self.CommutatorForInd(numPairs, ...
                                          ind - (numPairs - 1) * numPairs);
            end
        end
        
        function cdiisCoeffVec = InitialCoeffVec(~, tensorUse)
            numPairs = size(tensorUse, 1);
            onesVec = ones(numPairs, 1);
            hessian = zeros(numPairs, numPairs);
            for i = 1:numPairs
                for j = 1:numPairs
                    hessian(i, j) = tensorUse(i, i, j, j);
                end
            end
            hessian = [hessian, onesVec; onesVec', 0];
            cdiisCoeffVec = hessian \ [zeros(numPairs,1); 1];
            cdiisCoeffVec = cdiisCoeffVec(1:end-1);
        end
        
        function coeffVec = FindCoeffVec(self)
            for numPairs = length(self.pairs):-1:2
                tensorUse = self.tensor(end-numPairs+1:end, ...
                                        end-numPairs+1:end, ...
                                        end-numPairs+1:end, ...
                                        end-numPairs+1:end);
                iniCoeffVec = self.InitialCoeffVec(tensorUse);
                coeffVec = self.NewtonSolver(tensorUse, iniCoeffVec);
                if isnan(sum(coeffVec))
                    disp('NaN found in coefficients; reducing tensor size')
                else
                    return;
                end
            end
            coeffVec = 1;
        end
        
        function [grad, hess] = GradHess(~, tensorGrad, tensorHess, coeffVec)
            numPairs = length(coeffVec);
            
            grad = reshape(tensorGrad, [], numPairs) * coeffVec;
            grad = reshape(grad, [], numPairs) * coeffVec;
            grad = reshape(grad, [], numPairs) * coeffVec;
            
            hess = reshape(tensorHess, [], numPairs) * coeffVec;
            hess = reshape(hess, [], numPairs) * coeffVec;
            hess = reshape(hess, [], numPairs);
        end
        
        function coeffVec = NewtonSolver(self, tensorUse, coeffVec)
            tensorGrad = tensorUse + permute(tensorUse, [2 1 3 4]);
            tensorHess = tensorUse + permute(tensorUse, [1 3 2 4]) ...
                                   + permute(tensorUse, [1 4 2 3]) ...
                                   + permute(tensorUse, [2 1 3 4]) ...
                                   + permute(tensorUse, [2 3 1 4]) ...
                                   + permute(tensorUse, [2 4 1 3]);
            lambda = 0;
            for iter = 1:self.maxIterNewton
                [grad, hess] = self.GradHess(tensorGrad, tensorHess, coeffVec);
                gradL = [grad + lambda; 0];
                onesVec = ones(length(coeffVec), 1);
                hessL = [hess, onesVec; onesVec', 0];
                if rcond(hessL) < self.minHessRcond
                    disp('Inversion failed')
                    return;
                end
                coeffsAndLambda = [coeffVec; lambda];
                coeffsAndLambda = coeffsAndLambda - hessL \ gradL;
                coeffVec = coeffsAndLambda(1:end-1);
                lambda = coeffsAndLambda(end);
                if sqrt(mean(gradL.^2)) < self.gradNormThres
                    return;
                end
            end
            
            if iter > 199
                disp('Not converged');
            end
        end
        
    end
    
end
