function matrix = G09ReadMatrix(type)
if(strcmpi(type, 'ecpInt'))
    beginning = ' ECP Integrals: ';
    ending = ' SVDSVc ';
elseif(strcmpi(type, 'HarrisGuessMOAlpha'))
    beginning = ' Guess MO coefficients \(alpha\):';
    ending = ' Guess MO coefficients \(beta\):';
elseif(strcmpi(type, 'HarrisGuessMOBeta'))
    beginning = ' Guess MO coefficients \(beta\):';
    ending = ' Initial guess ';
end

matrix = [];
logFile = fopen('temp.log');
currLine = '';
while(ischar(currLine))
    blocks = {};
    iBlock = 0;
    if(~isempty(regexp(currLine, beginning, 'ONCE')))
        readLine = fgetl(logFile);
        while(isempty(regexp(readLine, ending, 'ONCE')))
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
        matrix = currMat;
        break;
    end
    currLine = fgetl(logFile);
end
fclose(logFile);

end
