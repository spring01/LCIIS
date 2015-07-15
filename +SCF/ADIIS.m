classdef ADIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        densVectors;
        
    end
    
    methods
        
        function obj = ADIIS(initVector, numVectors, type)
            if(nargin < 2)
                numVectors = 5;
            end
            if(nargin < 3)
                type = 'r';
            end
            lenVec = numel(initVector);
            if(strcmpi(type, 'r'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
            elseif(strcmpi(type, 'u'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
                obj.fockVectors{2} = zeros(lenVec, numVectors);
                obj.densVectors{2} = zeros(lenVec, numVectors);
            end
        end
        
        function Push(obj, newFockVector, newDensVector)
            for spin = 1:length(obj.fockVectors)
                % push new Fock in
                obj.fockVectors{spin}(:, 1:end-1) = obj.fockVectors{spin}(:, 2:end);
                obj.fockVectors{spin}(:, end) = newFockVector(:, spin);
                
                % push new density in
                obj.densVectors{spin}(:, 1:end-1) = obj.densVectors{spin}(:, 2:end);
                obj.densVectors{spin}(:, end) = newDensVector(:, spin);
            end
        end
        
        function [optFockVector, coeffs, useFockVectors] = OptFockVector(obj)
            useFockVectors = cell(1, length(obj.fockVectors));
            optFockVector = zeros(size(obj.fockVectors{1}, 1), length(obj.fockVectors));
            numVectors = sum(sum(obj.densVectors{1}.^2) ~= 0);
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
            
            errorFockVectors = cell(1, length(obj.fockVectors));
            errorDensVectors = cell(1, length(obj.fockVectors));
            % compute errors wrt. the latest Fock or density
            for spin = 1:length(obj.fockVectors)
                errorFockVectors{spin} = useFockVectors{spin} - ...
                    repmat(useFockVectors{spin}(:,end), 1, numVectors);
                errorDensVectors{spin} = useDensVectors{spin} - ...
                    repmat(useDensVectors{spin}(:,end), 1, numVectors);
            end
            
            % first order term and Hessian
            firstOrder = zeros(numVectors, 1);
            hessian = zeros(numVectors, numVectors);
            for spin = 1:length(obj.fockVectors)
                firstOrder = firstOrder + 2.*(useFockVectors{spin}(:, end)'*errorDensVectors{spin})';
                hessian = hessian + errorFockVectors{spin}'*errorDensVectors{spin};
            end
            hessian = hessian + hessian'; % multiply Hessian by 2 and cancels numerical error
            
            options = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
            coeffs = quadprog(real(hessian), real(firstOrder), ...
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