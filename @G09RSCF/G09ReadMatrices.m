function matrices = G09ReadMatrices(types)
beginning = cell(1, length(types));
ending = cell(1, length(types));
for iMat = 1:length(types)
    if(strcmpi(types{iMat}, 'overlap'))
        beginning{iMat} = ' *** Overlap *** ';
        ending{iMat} = ' *** Kinetic Energy *** ';
    elseif(strcmpi(types{iMat}, 'coreHamilt'))
%         beginning{iMat} = ' ****** Core Hamiltonian ****** ';
%         ending{iMat} = ' SVDSVc ';
%         ending{iMat} = ' Symmetry operations used ';
        beginning{iMat} = ' Core hamiltonian \(alpha\):';
        ending{iMat} = ' Core hamiltonian \(beta\):';
    elseif(strcmpi(types{iMat}, 'fockAlpha'))
        beginning{iMat} = ' Fock matrix \(alpha\):';
        ending{iMat} = ' Fock matrix \(beta\):';
    elseif(strcmpi(types{iMat}, 'fockBeta'))
        beginning{iMat} = ' Fock matrix \(beta\):';
        ending{iMat} = ' E= ';
    elseif(strcmpi(types{iMat}, 'ecpInt'))
        beginning{iMat} = ' ECP Integrals: ';
        ending{iMat} = ' SVDSVc ';
    elseif(strcmpi(types{iMat}, 'orbAlpha'))
        beginning{iMat} = ' Alpha MO coefficients ';
        ending{iMat} = ' Alpha density matrix ';
    elseif(strcmpi(types{iMat}, 'orbBeta'))
        beginning{iMat} = ' Beta MO coefficients ';
        ending{iMat} = ' Beta density matrix ';
    elseif(strcmpi(types{iMat}, 'HarrisGuessMO'))
        beginning{iMat} = ' Guess MO coefficients:';
        ending{iMat} = ' Leave Link ';
    elseif(strcmpi(types{iMat}, 'HarrisGuessMOAlpha'))
        beginning{iMat} = ' Guess MO coefficients \(alpha\):';
        ending{iMat} = ' Guess MO coefficients \(beta\):';
    elseif(strcmpi(types{iMat}, 'HarrisGuessMOBeta'))
        beginning{iMat} = ' Guess MO coefficients \(beta\):';
        ending{iMat} = ' Initial guess ';
    end
end

matrices = {};

logFile = fopen('temp.log');
currLine = '';
while(ischar(currLine))
    for iMat = 1:length(types)
        blocks = {};
        iBlock = 0;
        if(~isempty(regexp(currLine, beginning{iMat}, 'ONCE')))
            readLine = fgetl(logFile);
            while(isempty(regexp(readLine, ending{iMat}, 'ONCE')))
                allMatches = regexp(readLine, '[+-]?[0-9]+.[0-9]+D[+-][0-9][0-9]', 'match');
                if(~isempty(allMatches)) % rowNum with doubles
                    numsInALine = zeros(1, numOfNums);
                    for iMatch = 1:length(allMatches)
                        numsInALine(iMatch) = str2num(allMatches{iMatch}); %#ok
                    end
                    rowNum = str2double(regexp(readLine, '[0-9]+', 'match', 'ONCE'));
                    currentBlock(rowNum, :) = numsInALine; %#ok
                    blocks{iBlock} = currentBlock; %#ok
                else % integers
                    intMatches = regexp(readLine, '[0-9]+', 'match');
                    numOfNums = length(intMatches);
                    iBlock = iBlock + 1;
                    currentBlock = zeros(0, numOfNums);
                end
                readLine = fgetl(logFile);
            end
            currMat = [blocks{:}];
            if(triu(currMat, 1) == 0)
                currMat = currMat + currMat' - diag(diag(currMat));
            end
            matrices{iMat} = currMat; %#ok
            
            currLine = readLine;
        end
    end
    currLine = fgetl(logFile);
end
fclose(logFile);

end
