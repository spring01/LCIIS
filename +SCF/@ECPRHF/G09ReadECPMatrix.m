function ecpMat = G09ReadECPMatrix()
logFile = fopen('temp.log');
currLine = '';
ecpMat = [];
while(ischar(currLine))
    blocks = {};
    iBlock = 0;
    if(~isempty(regexp(currLine, ' ECP Integrals: ', 'ONCE')))
        readLine = fgetl(logFile);
        while(isempty(regexp(readLine, ' SVDSVc ', 'ONCE')))
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
        ecpMat = currMat;
        break;
    end
    currLine = fgetl(logFile);
end
fclose(logFile);

end
