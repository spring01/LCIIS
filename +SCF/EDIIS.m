classdef EDIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        densVectors;
        energies;
        
    end
    
    methods
        
        function obj = EDIIS(oeiVector, numVectors, type)
            if(nargin < 2)
                numVectors = 5;
            end
            if(nargin < 3)
                type = 'r';
            end
            oeiVector = reshape(oeiVector, [], 1);
            lenVec = length(oeiVector);
            if(strcmpi(type, 'r'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
            elseif(strcmpi(type, 'u'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
                obj.fockVectors{2} = zeros(lenVec, numVectors);
                obj.densVectors{2} = zeros(lenVec, numVectors);
            end
            obj.energies = zeros(numVectors, 1);
        end
        
        function Push(obj, newFockVector, newDensVector, newElecEnergy)
            % push new Fock and density in
            for spin = 1:length(obj.fockVectors)
                obj.fockVectors{spin}(:, 1:end-1) = obj.fockVectors{spin}(:, 2:end);
                obj.fockVectors{spin}(:, end) = newFockVector(:, spin);
                
                obj.densVectors{spin}(:, 1:end-1) = obj.densVectors{spin}(:, 2:end);
                obj.densVectors{spin}(:, end) = newDensVector(:, spin);
            end
            
            % push new energy in
            obj.energies(1:end-1) = obj.energies(2:end);
            obj.energies(end) = newElecEnergy;
        end
        
        function [optFockVector, coeffs, useFockVectors] = OptFockVector(obj)
            numVectors = sum(sum(obj.densVectors{1}.^2) ~= 0);
            [optFockVector, coeffs, useFockVectors] = obj.SolveForNumVectors(numVectors);
        end
        
        function [optFockVector, coeffs, useFockVectors] = SolveForNumVectors(obj, numVectors)
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
            for spin = 1:length(obj.fockVectors)
                useFockVectors{spin} = obj.fockVectors{spin}(:, end-numVectors+1:end);
                useDensVectors{spin} = obj.densVectors{spin}(:, end-numVectors+1:end);
            end
            energies_ = obj.energies(end-numVectors+1:end);
            
            hessian = zeros(numVectors, numVectors);
            for i = 1:numVectors
                for j = 1:numVectors
                    errorFock = [];
                    errorDens = [];
                    for spin = 1:length(obj.fockVectors)
                        errorFock = [errorFock; useFockVectors{spin}(:, i) - useFockVectors{spin}(:, j)]; %#ok
                        errorDens = [errorDens; useDensVectors{spin}(:, i) - useDensVectors{spin}(:, j)]; %#ok
                    end
                    hessian(j, i) = errorFock' * errorDens;
                end
            end
            hessian = -hessian ./ length(obj.fockVectors);
            
            options = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
            coeffs = quadprog(hessian, energies_, ...
                -eye(numVectors), zeros(numVectors, 1), ...
                ones(1,numVectors), 1, ...
                [], [], ...
                [], options);
            
            for spin = 1:length(obj.fockVectors)
                optFockVector(:, spin) = useFockVectors{spin} * coeffs;
            end
            
            if(length(useFockVectors) == 1)
                useFockVectors = useFockVectors{1};
            end
            
%             disp(coeffs');
        end
        
    end
    
end