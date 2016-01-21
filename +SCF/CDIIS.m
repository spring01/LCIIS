classdef CDIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        errorVectors;
        
        S_Half;
        inv_S_Half;
        
        minHessRcond = 1e-20;
        
    end
    
    methods
        
        function obj = CDIIS(overlapMatrix, numVectors, type)
            if nargin < 2
                numVectors = 5;
            end
            if nargin < 3
                type = 'r';
            end
            lenVec = numel(overlapMatrix);
            if strcmpi(type, 'r')
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.errorVectors = zeros(lenVec, numVectors);
            elseif(strcmpi(type, 'u'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.fockVectors{2} = zeros(lenVec, numVectors);
                obj.errorVectors = zeros(2*lenVec, numVectors);
            end
            obj.S_Half = sqrtm(overlapMatrix);
            obj.inv_S_Half = inv(obj.S_Half);
        end
        
        function Push(obj, newFockVector, newDensVector)
            % push new Fock in
            for spin = 1:length(obj.fockVectors)
                obj.fockVectors{spin}(:, 1:end-1) = obj.fockVectors{spin}(:, 2:end);
                obj.fockVectors{spin}(:, end) = newFockVector(:, spin);
            end
            
            % push new commutator error in
            obj.errorVectors(:, 1:end-1) = obj.errorVectors(:, 2:end);
            errorVec = [];
            for spin = 1:length(obj.fockVectors)
                FtDt = obj.inv_S_Half ...
                    * reshape(newFockVector(:, spin), sqrt(length(newFockVector(:, spin))), []) ...
                    * reshape(newDensVector(:, spin), sqrt(length(newDensVector(:, spin))), []) ...
                    * obj.S_Half;
                errorVec = [errorVec; reshape(FtDt - FtDt', [], 1)]; %#ok
            end
            obj.errorVectors(:, end) = errorVec;
        end
        
        function [optFockVector, coeffs, useFockVectors] = OptFockVector(obj, numVectors)
            if nargin < 2
                numVectors = sum(sum(obj.errorVectors.^2) ~= 0);
            end
            useFockVectors = cell(1, length(obj.fockVectors));
            optFockVector = zeros(size(obj.fockVectors{1}, 1), length(obj.fockVectors));
            if numVectors == 0 || numVectors == 1
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
            useErrorVectors = obj.errorVectors(:, end-numVectors+1:end);
            for spin = 1:length(obj.fockVectors)
                useFockVectors{spin} = obj.fockVectors{spin}(:, end-numVectors+1:end);
            end
            
            onesVec = ones(numVectors, 1);
            hessian = [ ...
                useErrorVectors'*useErrorVectors, onesVec; ...
                onesVec', 0];
            hessian = (hessian + hessian') ./ 2;
            if rcond(hessian) < obj.minHessRcond || isnan(rcond(hessian))
                disp('Inversion failed')
                diisCoefficients = NaN .* [zeros(numVectors,1); 1];
            else
                diisCoefficients = linsolve(hessian, [zeros(numVectors,1); 1]);
            end
            coeffs = diisCoefficients(1:end-1);
            
            for spin = 1:length(obj.fockVectors)
                optFockVector(:, spin) = useFockVectors{spin} * coeffs;
            end
            
            if length(useFockVectors) == 1
                useFockVectors = useFockVectors{1};
            end
            
            if isnan(coeffs)
                numVectors = numVectors - 1;
                [optFockVector, coeffs, useFockVectors] = obj.OptFockVector(numVectors);
            end
        end
        
        function optFockVector = CalcFockVec(obj, coeffs, useFockVectors)
            optFockVector = zeros(size(obj.fockVectors{1}, 1), length(obj.fockVectors));
            if length(obj.fockVectors) == 1
                newCell{1} = useFockVectors;
                useFockVectors = newCell;
            end
            for spin = 1:length(obj.fockVectors)
                optFockVector(:, spin) = useFockVectors{spin} * coeffs;
            end
        end
        
        function maxError = MaxError(obj)
            maxError = max(abs(obj.errorVectors(:, end)));
        end
        
    end
    
end