function orb = G09OrbToPsi4Orb(self, orb)
shellNfuncs = self.matpsi.BasisSet_ShellNumFunctions();
shellStartFuncs = cumsum([1 shellNfuncs]);
shellStartFuncs = shellStartFuncs(1:end-1);
basisSetIsSpherical = self.matpsi.BasisSet_IsSpherical();
for i = 1:length(shellNfuncs)
    offset = shellStartFuncs(i) - 1;
    renormInd = [];
    renormBy = 1;
    srcInd = [];
    tgtInd = [];
    if shellNfuncs(i) == 3 && basisSetIsSpherical
        % 3p spherical; need to change [x y z] -> [z x y]
        srcInd = [3, 1, 2] + offset;
        tgtInd = (1:3) + offset;
    elseif shellNfuncs(i) == 6 && ~basisSetIsSpherical
        % 6d cartesian; [xx yy zz xy xz yz] -> [xx xy xz yy yz zz]
        renormInd = (4:6) + offset;
        renormBy = sqrt(3);
        srcInd = [1, 4, 5, 2, 6, 3] + offset;
        tgtInd = (1:6) + offset;
    elseif shellNfuncs(i) == 10 && ~basisSetIsSpherical
        % 10f cartesian
        % [xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz] ->
        % [xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz]
        renormInd = (4:10) + offset;
        renormBy = repmat([sqrt(5) * ones(1, 6), sqrt(15)]', ...
            1, size(orb, 2));
        srcInd = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3] + offset;
        tgtInd = (1:10) + offset;
    elseif shellNfuncs(i) == 15 && ~basisSetIsSpherical
        % 15g cartesian; reverse order
        renormInd = (2:14) + offset;
        renormBy = repmat([sqrt(7), sqrt(35/3), sqrt(7), ...
            1, sqrt(7), sqrt(35), sqrt(35), ...
            sqrt(7), sqrt(35/3), sqrt(35), ...
            sqrt(35/3), sqrt(7), sqrt(7)]', ...
            1, size(orb, 2));
        srcInd = (15:-1:1) + offset;
        tgtInd = (1:15) + offset;
    end
    orb(renormInd, :) = orb(renormInd, :) .* renormBy;
    orb(tgtInd, :) = orb(srcInd, :);
end
end