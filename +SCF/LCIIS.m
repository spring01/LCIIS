classdef LCIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        densVectors;
        S_Half;
        inv_S_Half;
        
        tempTensorForGrad;
        tempTensorForHess;
        
        minHessRcond = 1e-25;
        gradNormThres = 1e-12;
        
    end
    
    methods
        
        function obj = LCIIS(overlapMatrix, numVectors, type)
            if(nargin < 2)
                numVectors = 5;
            end
            if(nargin < 3)
                type = 'r';
            end
            lenVec = numel(overlapMatrix);
            if(strcmpi(type, 'r'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
            elseif(strcmpi(type, 'u'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
                obj.fockVectors{2} = zeros(lenVec, numVectors);
                obj.densVectors{2} = zeros(lenVec, numVectors);
            end
            obj.S_Half = sqrtm(overlapMatrix);
            obj.inv_S_Half = inv(obj.S_Half);
        end
        
        function Push(obj, newFockVector, newDensVector)
            for spin = 1:length(obj.fockVectors)
                % push Fock
                obj.fockVectors{spin}(:, 1:end-1) = obj.fockVectors{spin}(:, 2:end);
                obj.fockVectors{spin}(:, end) = newFockVector(:, spin);
                
                % push density
                obj.densVectors{spin}(:, 1:end-1) = obj.densVectors{spin}(:, 2:end);
                obj.densVectors{spin}(:, end) = newDensVector(:, spin);
            end
        end
        
        function [optFockVector, coeffs, useFockVectors] = OptFockVector(obj, numVectors)
            if(nargin < 2)
                numVectors = sum(sum(obj.densVectors{1}.^2) ~= 0);
            end
            useFockVectors = cell(1, length(obj.fockVectors));
            optFockVector = zeros(size(obj.fockVectors{1}, 1), length(obj.fockVectors));
            if(numVectors == 0 || numVectors == 1)
                for spin = 1:length(obj.fockVectors)
                    optFockVector(:, spin) = obj.fockVectors{spin}(:, end);
                    useFockVectors{spin} = obj.fockVectors{spin}(:, end);
                end
                if(length(useFockVectors) == 1)
                    useFockVectors = useFockVectors{1};
                end
                coeffs = 1;
                return;
            end
            useDensVectors = cell(1, length(obj.fockVectors));
            tensor = zeros(numVectors, numVectors, numVectors, numVectors);
            for spin = 1:length(obj.fockVectors)
                useFockVectors{spin} = obj.fockVectors{spin}(:, end-numVectors+1:end);
                useDensVectors{spin} = obj.densVectors{spin}(:, end-numVectors+1:end);
                tensor = tensor + obj.CalcTensor(useFockVectors{spin}, useDensVectors{spin});
            end
            
            % start from CDIIS coefficients
            onesVec = ones(numVectors, 1);
            hessian = zeros(numVectors, numVectors);
            for errorI = 1:numVectors
                for errorJ = 1:numVectors
                    hessian(errorI, errorJ) = tensor(errorI, errorI, errorJ, errorJ);
                end
            end
            hessian = [hessian, onesVec; onesVec', 0];
            diisCoefficients = hessian \ [zeros(numVectors,1); 1];
            iniCoeffs = diisCoefficients(1:end-1);
            
            coeffs = obj.NewtonSolver(tensor, iniCoeffs);
            
            for spin = 1:length(obj.fockVectors)
                optFockVector(:, spin) = useFockVectors{spin} * coeffs;
            end
            
            if(length(useFockVectors) == 1)
                useFockVectors = useFockVectors{1};
            end
            
            if(isnan(coeffs))
                numVectors = numVectors - 1;
                [optFockVector, coeffs, useFockVectors] = obj.OptFockVector(numVectors);
            end
        end
                
    end
    
    methods (Access = private)
        
        function tensor = CalcTensor(obj, useFockVectors, useDensVectors)
            nbf = sqrt(size(useFockVectors, 1));
            numVectors = size(useFockVectors, 2);
            error = zeros(numVectors^2, nbf^2);
            for dummyI = 1:numVectors
                fockI = reshape(useFockVectors(:, dummyI), nbf, []);
                for dummyJ = 1:numVectors
                    densJ = reshape(useDensVectors(:, dummyJ), nbf, []);
                    FtDt = obj.inv_S_Half * fockI * densJ * obj.S_Half;
                    error((dummyI-1)*numVectors+dummyJ, :) = reshape(FtDt - FtDt', [], 1);
                end
            end
            tensor = reshape(error*error', numVectors, numVectors, numVectors, numVectors);
        end
        
        function [grad, hess] = GradHess(obj, coeffs)
            nVecs = length(coeffs);
            
            grad = reshape(obj.tempTensorForGrad, [], nVecs) * coeffs;
            grad = reshape(grad, [], nVecs) * coeffs;
            grad = reshape(grad, [], nVecs) * coeffs;
            
            hess = reshape(obj.tempTensorForHess, [], nVecs) * coeffs;
            hess = reshape(hess, [], nVecs) * coeffs;
            hess = reshape(hess, [], nVecs);
        end
        
        function coeffs = NewtonSolver(obj, tensor, iniCoeffs)
            coeffs = iniCoeffs;
            obj.tempTensorForGrad = 2*(tensor + permute(tensor, [2 1 3 4]));
%             obj.tempTensorForHess = ...
%                 2*(tensor + permute(tensor, [1 3 2 4]) + permute(tensor, [1 4 2 3]) + ...
%                 permute(tensor, [2 1 3 4]) + permute(tensor, [2 3 1 4]) + permute(tensor, [2 4 1 3]));
            obj.tempTensorForHess = ...
                2*(tensor + ipermute(tensor, [1 3 2 4]) + ipermute(tensor, [1 3 4 2]) + ...
                ipermute(tensor, [2 1 3 4]) + ipermute(tensor, [3 1 2 4]) + ipermute(tensor, [3 1 4 2]));
            lambda = 0;
            for iter = 1:200
                [grad, hess] = obj.GradHess(coeffs);
                gradL = [grad + lambda; 0];
                onesVec = ones(length(coeffs), 1);
                hessL = [hess, onesVec; onesVec', 0];
                if(rcond(hessL) < obj.minHessRcond)
                    disp('Inversion failed')
                    return;
                end
                coeffsAndLambda = [coeffs; lambda];
                coeffsAndLambda = coeffsAndLambda - hessL \ gradL;
                coeffs = coeffsAndLambda(1:end-1);
                lambda = coeffsAndLambda(end);
                if(sqrt(mean(gradL.^2)) < obj.gradNormThres)
                    return;
                end
            end
            
            if(iter > 199)
                disp('Not converged');
            end
        end
        
    end
    
end
