classdef RLCIIS < handle
    
    properties (Access = protected)
        
        pairs;
        comms;
        
    end
    
    properties (Access = private)
        
        sqrtmS;
        invSqrtmS;
        
        maxNumPairs;
        firstUse;
        pairsFull;
        bigMat;
        
        maxIterNewton = 200;
        minHessRcond = 1e-20;
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
            self.firstUse = true;
        end
        
        function Push(self, fockVector, densVector)
            newPair = self.CreatePair(fockVector, densVector);
            if length(self.pairs) == self.maxNumPairs
                self.pairs(1:end-1) = self.pairs(2:end);
                self.pairs{end} = newPair;
                self.pairsFull = true;
            else
                self.pairs{end+1} = newPair;
                self.pairsFull = false;
            end
            
            if ~self.firstUse
                if self.pairsFull
                    self.PreUpdateFull();
                else
                    self.PreUpdateNotFull();
                end
                self.UpdateComms();
                self.UpdateBigMat();
            end
        end
        
        function [optFockVector, coeffVec, useFockVectors] = OptFockVector(self)
            if self.firstUse
                self.CalcTensor();
                self.firstUse = false;
            end
            coeffVec = self.FindCoeffVec();
            [optFockVector, useFockVectors] = self.ExtrapolateFock(coeffVec);
        end
        
    end
    
    methods (Access = protected)
        
        function pair = CreatePair(self, fockVector, densVector)
            % This function defines interface.
            % In restricted SCF, newFockVector stores the Fock matrix, and
            % newDensVector stores the density matrix, both as an nbf^2 by 1
            % vector.
            pair.fockVector = fockVector;
            pair.fockOrthoMatrix = self.FockOrthoMat(fockVector);
            pair.densOrthoMatrix = self.DensOrthoMat(densVector);
        end
        
        function fockOrthoMat = FockOrthoMat(self, singleFockVector)
            nbf = size(self.invSqrtmS, 1);
            fockOrthoMat = self.invSqrtmS * reshape(singleFockVector, nbf, nbf);
        end
        
        function densOrthoMat = DensOrthoMat(self, singleDensVector)
            nbf = size(self.sqrtmS, 1);
            densOrthoMat = reshape(singleDensVector, nbf, nbf) * self.sqrtmS;
        end
        
        function comm = CommBetween(self, ind1, ind2)
            comm = self.Commutator(self.pairs{ind1}.fockOrthoMatrix, ...
                                   self.pairs{ind2}.densOrthoMatrix);
        end
        
        function comm = Commutator(~, fockOrthoMatrix, densOrthoMatrix)
            fockOrthoDensOrtho = fockOrthoMatrix * densOrthoMatrix;
            comm = reshape(fockOrthoDensOrtho - fockOrthoDensOrtho', [], 1);
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
            for i = 1:numPairs
                for j = 1:numPairs
                    commInd = (i - 1) * numPairs + j;
                    self.comms(:, commInd) = self.CommBetween(i, j);
                end
            end
            self.bigMat = self.comms' * self.comms;
        end
        
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
        
        function cdiisCoeffVec = InitialCoeffVec(self, tensorUse)
            numPairs = size(tensorUse, 1);
            onesVec = ones(numPairs, 1);
            hess = zeros(numPairs, numPairs);
            for i = 1:numPairs
                for j = 1:numPairs
                    hess(i, j) = tensorUse(i, i, j, j);
                end
            end
            hess = (hess + hess') ./ 2;
            hessL = [hess, onesVec; onesVec', 0];
            rcondHessL = rcond(hessL);
            if rcondHessL < self.minHessRcond || isnan(rcondHessL)
                disp('Inversion failed');
                cdiisCoeffVec = NaN .* ones(numPairs, 1);
                return;
            end
            cdiisCoeffVec = linsolve(hessL, [zeros(numPairs, 1); 1]);
            cdiisCoeffVec = cdiisCoeffVec(1:(end - 1));
        end
        
        function coeffVec = FindCoeffVec(self)
            numPairs = length(self.pairs);
            allInd = 1:numPairs^2;
            tensor = reshape(self.bigMat(allInd, allInd), ...
                             numPairs, numPairs, numPairs, numPairs);
            for useNumPairs = 1:(numPairs - 1)
                useInd = useNumPairs:numPairs;
                tensorUse = tensor(useInd, useInd, useInd, useInd);
                iniCoeffVec = self.InitialCoeffVec(tensorUse);
                if isnan(sum(iniCoeffVec))
                    disp('NaN found in iniCoeffVec; reducing tensor size');
                    continue;
                end
                coeffVec = self.NewtonSolver(tensorUse, iniCoeffVec);
                if isnan(sum(coeffVec))
                    disp('NaN found in coefficients; reducing tensor size');
                    continue;
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
        
        function coeffVec = NewtonSolver(self, tensorUse, iniCoeffVec)
            coeffVec = iniCoeffVec;
            tensorGrad = tensorUse + permute(tensorUse, [2 1 3 4]);
            tensorHess = tensorUse + permute(tensorUse, [1 3 2 4]) ...
                                   + permute(tensorUse, [1 4 2 3]) ...
                                   + permute(tensorUse, [2 1 3 4]) ...
                                   + permute(tensorUse, [2 3 1 4]) ...
                                   + permute(tensorUse, [2 4 1 3]);
            lambda = 0;
            for iter = 1:self.maxIterNewton
                oldCoeffVec = coeffVec;
                [grad, hess] = self.GradHess(tensorGrad, tensorHess, coeffVec);
                hess = (hess + hess') ./ 2;
                gradL = [grad + lambda; 0];
                onesVec = ones(length(coeffVec), 1);
                hessL = [hess, onesVec; onesVec', 0];
                rcondHessL = rcond(hessL);
                if rcondHessL < self.minHessRcond || isnan(rcondHessL)
                    disp('Inversion failed')
                    coeffVec = oldCoeffVec;
                    return;
                end
                coeffsAndLambda = [coeffVec; lambda];
                coeffsAndLambda = coeffsAndLambda - linsolve(hessL, gradL);
                coeffVec = coeffsAndLambda(1:end-1);
                lambda = coeffsAndLambda(end);
                if sqrt(mean(gradL.^2)) < self.gradNormThres
                    return;
                end
            end
            
            if iter > self.maxIterNewton - 1
                disp('Not converged');
            end
        end
        
    end
    
end
