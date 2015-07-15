classdef Molecule < handle
    
    properties (SetAccess = private)
        
        cartesian;
        charge = 0;
        multiplicity = 1;
        
    end
    
    
    properties (Constant, Hidden = true)
        
        Bohr2Angstrom = 0.529177249;
        
    end
        
    methods
        
        function obj = Molecule(cartesian, charge, multiplicity)
            if(nargin > 1)
                obj.charge = charge;
            end
            if(nargin > 2)
                obj.multiplicity = multiplicity;
            end
            obj.cartesian = cartesian;
        end
        
        function numAtoms = NumAtoms(obj)
            numAtoms = size(obj.cartesian, 1);
        end
        
        function numElectrons = NumElectrons(obj)
            numElectrons = sum(obj.cartesian(:, 1)) - obj.charge;
        end
        
        function Plot(obj)
            tempFileName = [tempname(), '.xyz'];
            tempFile = fopen(tempFileName, 'w');
            xyzStr = [num2str(size(obj.cartesian,1)), char(10), 'temp', char(10), obj.MoleculeString()];
            fprintf(tempFile, xyzStr);
            system(['jmol ', tempFileName]);
            system(['rm ', tempFileName]);
        end
        
    end
    
    methods (Access = private)
        
        function molString = MoleculeString(obj)
            numAtoms = size(obj.cartesian,1);
            symbols = repmat(' ', numAtoms, 5);
            periodicTable = obj.PeriodicTable();
            for iatom = 1:numAtoms
                atomSymbol = periodicTable{obj.cartesian(iatom, 1)};
                symbols(iatom, 1:length(atomSymbol)) = atomSymbol;
            end
            formSpec = '%.10f';
            spaces = repmat('  ', numAtoms, 1);
            molString = [symbols, num2str(obj.cartesian(:, 2), formSpec), ...
                spaces, num2str(obj.cartesian(:, 3), formSpec), ...
                spaces, num2str(obj.cartesian(:, 4), formSpec), ...
                repmat(char(10), numAtoms, 1)];
            molString = reshape(molString', 1, []);
        end
        
        function periodicTable = PeriodicTable(~)
            periodicTable = { ...
                'H',                                                                                                  'He', ...
                'Li', 'Be',                                                             'B' , 'C' , 'N' , 'O' , 'F' , 'Ne', ...
                'Na', 'Mg',                                                             'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', ...
                'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', ...
                'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe', ...
                'Cs', 'Ba',  ...
                            'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', ...
                                  'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', ...
                'Fr', 'Ra', ...
                            'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', ...
                                  'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn','Uut', 'Fl','Uup', 'Lv','Uus','Uuo'};
        end
        
    end
    
end
