classdef RLCIIS < handle
    
    properties (Access = protected)
        
        pairs;
        comms;
        
    end
    
    properties (Access = private)
        
        sqrtmS;
        invSqrtmS;
        maxNumPairs;
        pairsFull;
        bigMat;
        
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
            self.comms = zeros(numel(overlapMatrix), maxNumPairs^2);
            self.bigMat = zeros(maxNumPairs^2, maxNumPairs^2);
            self.sqrtmS = sqrtm(overlapMatrix);
            self.invSqrtmS = inv(self.sqrtmS);
        end
        
        function Push(self, newFockVector, newDensVector)
            newPair.fockVector = newFockVector;
            newPair.densVector = newDensVector;
            if length(self.pairs) == size(self.comms, 1)
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
                self.PreUpdateFull();
            else
                self.PreUpdateNotFull();
            end
            self.UpdateComms();
            self.UpdateBigMat();
            coeffVec = self.FindCoeffVec();
            [optFockVector, useFockVectors] = self.ExtrapolateFock(coeffVec);
        end
        
    end
    
    methods (Access = protected)
        
        function comm = CommBetween(self, ind1, ind2)
            comm = self.Commutator(self.pairs{ind1}.fockVector, ...
                                   self.pairs{ind2}.densVector);
        end
        
        function comm = Commutator(self, singleFockVector, singleDensVector)
            nbf = size(self.sqrtmS, 1);
            fockMat = reshape(singleFockVector, nbf, nbf);
            densMat = reshape(singleDensVector, nbf, nbf);
            fockDens = self.invSqrtmS * fockMat * densMat * self.sqrtmS;
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
        
        function PreUpdateFull(self)
            numPairs = self.maxNumPairs;
            for i = 2:numPairs
                srcInd = ((i - 1) * numPairs + 2):(i * numPairs);
                tgtInd = srcInd - numPairs - 1;
                self.comms(:, tgtInd) = self.comms(:, srcInd);
                self.bigMat(tgtInd, :) = self.bigMat(srcInd, :);
                self.bigMat(:, tgtInd) = self.bigMat(:, srcInd);
            end
        end
        
        function PreUpdateNotFull(self)
            numPairs = length(self.pairs);
            for i = (numPairs - 1):-1:2
                srcInd = ((i - 1) * (numPairs - 1) + 1):(i * (numPairs - 1));
                tgtInd = srcInd + i - 1;
                self.comms(:, tgtInd) = self.comms(:, srcInd);
                self.bigMat(tgtInd, :) = self.bigMat(srcInd, :);
                self.bigMat(:, tgtInd) = self.bigMat(:, srcInd);
            end
        end
        
        function UpdateComms(self)
            numPairs = length(self.pairs);
            for indPair = 1:(numPairs-1)
                indComm = indPair * numPairs;
                self.comms(:, indComm) = self.CommBetween(indPair, numPairs);
            end
            
            for indComm = ((numPairs - 1) * numPairs + 1):numPairs^2
                indPair = indComm - (numPairs - 1) * numPairs;
                self.comms(:, indComm) = self.CommBetween(numPairs, indPair);
            end
        end
        
        function UpdateBigMat(self)
            numPairs = length(self.pairs);
            updInd = [numPairs * (1:(numPairs - 1)), ...
                    ((numPairs - 1) * numPairs + 1):numPairs^2];
            allInd = 1:numPairs^2;
            com = self.comms;
            self.bigMat(updInd, allInd) = com(:, updInd)' * com(:, allInd);
            self.bigMat(allInd, updInd) = self.bigMat(updInd, allInd)';
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
            numPairs = length(self.pairs);
            allInd = 1:numPairs^2;
            tensor = reshape(self.bigMat(allInd, allInd), numPairs*ones(1, 4));
            for useNumPairs = numPairs:-1:2
                useInd = numPairs-useNumPairs+1:numPairs;
                tensorUse = tensor(useInd, useInd, useInd, useInd);
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
