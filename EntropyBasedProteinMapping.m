%Read the protein fasta file
proteinSample = fastaread('.....\VP24.fasta');
%Conversion of struct to cell
proteinSample = struct2cell(proteinSample);
%Conversion of cell to mat (According to the protein sequence in the 1st
%place)
proteinSample = cell2mat(proteinSample(2));
%Determining the length of the protein sequence
maxbase = length(proteinSample)
%Creating an empty array
Z = []

%Conversion of the protein sequence with the proposed entropy-based method
for i=1:maxbase
    switch proteinSample(i)
        case('A')
            z = (seqwordcount(proteinSample,'A')/maxbase) * log(seqwordcount(proteinSample,'A')/maxbase);
            z = z *(seqwordcount(proteinSample,'A')/maxbase)^(1/log(seqwordcount(proteinSample,'A')/maxbase));
            Z=[Z z];
        case('R')
            z = (seqwordcount(proteinSample,'R')/maxbase) * log(seqwordcount(proteinSample,'R')/maxbase);
            z = z *(seqwordcount(proteinSample,'R')/maxbase)^(1/log(seqwordcount(proteinSample,'R')/maxbase));
            Z=[Z z]; 
        case('N')
            z = (seqwordcount(proteinSample,'N')/maxbase) * log(seqwordcount(proteinSample,'N')/maxbase);
            z = z *(seqwordcount(proteinSample,'N')/maxbase)^(1/log(seqwordcount(proteinSample,'N')/maxbase));
            Z=[Z z]; 
        case('D')
            z = (seqwordcount(proteinSample,'D')/maxbase) * log(seqwordcount(proteinSample,'D')/maxbase);
            z = z *(seqwordcount(proteinSample,'D')/maxbase)^(1/log(seqwordcount(proteinSample,'D')/maxbase));
            Z=[Z z]; 
        case('C')
            z = (seqwordcount(proteinSample,'C')/maxbase) * log(seqwordcount(proteinSample,'C')/maxbase);
            z = z *(seqwordcount(proteinSample,'C')/maxbase)^(1/log(seqwordcount(proteinSample,'C')/maxbase));
            Z=[Z z]; 
        case('Q')
            z = (seqwordcount(proteinSample,'Q')/maxbase) * log(seqwordcount(proteinSample,'Q')/maxbase);
            z = z *(seqwordcount(proteinSample,'Q')/maxbase)^(1/log(seqwordcount(proteinSample,'Q')/maxbase));
            Z=[Z z]; 
        case('E')
            z = (seqwordcount(proteinSample,'E')/maxbase) * log(seqwordcount(proteinSample,'E')/maxbase);
            z = z *(seqwordcount(proteinSample,'E')/maxbase)^(1/log(seqwordcount(proteinSample,'E')/maxbase));
            Z=[Z z]; 
        case('G')
            z = (seqwordcount(proteinSample,'G')/maxbase) * log(seqwordcount(proteinSample,'G')/maxbase);
            z = z *(seqwordcount(proteinSample,'G')/maxbase)^(1/log(seqwordcount(proteinSample,'G')/maxbase));
            Z=[Z z]; 
        case('H')
            z = (seqwordcount(proteinSample,'H')/maxbase) * log(seqwordcount(proteinSample,'H')/maxbase);
            z = z *(seqwordcount(proteinSample,'H')/maxbase)^(1/log(seqwordcount(proteinSample,'H')/maxbase));
            Z=[Z z]; 
        case('I')
            z = (seqwordcount(proteinSample,'I')/maxbase) * log(seqwordcount(proteinSample,'I')/maxbase);
            z = z *(seqwordcount(proteinSample,'I')/maxbase)^(1/log(seqwordcount(proteinSample,'I')/maxbase));
            Z=[Z z]; 
        case('L')
            z = (seqwordcount(proteinSample,'L')/maxbase) * log(seqwordcount(proteinSample,'L')/maxbase);
            z = z *(seqwordcount(proteinSample,'L')/maxbase)^(1/log(seqwordcount(proteinSample,'L')/maxbase));
            Z=[Z z]; 
        case('K')
            z = (seqwordcount(proteinSample,'K')/maxbase) * log(seqwordcount(proteinSample,'K')/maxbase);
            z = z *(seqwordcount(proteinSample,'K')/maxbase)^(1/log(seqwordcount(proteinSample,'K')/maxbase));
            Z=[Z z]; 
        case('M')
            z = (seqwordcount(proteinSample,'M')/maxbase) * log(seqwordcount(proteinSample,'M')/maxbase);
            z = z *(seqwordcount(proteinSample,'M')/maxbase)^(1/log(seqwordcount(proteinSample,'M')/maxbase));
            Z=[Z z]; 
        case('F')
            z = (seqwordcount(proteinSample,'F')/maxbase) * log(seqwordcount(proteinSample,'F')/maxbase);
            z = z *(seqwordcount(proteinSample,'F')/maxbase)^(1/log(seqwordcount(proteinSample,'F')/maxbase));
            Z=[Z z]; 
        case('P')
            z = (seqwordcount(proteinSample,'P')/maxbase) * log(seqwordcount(proteinSample,'P')/maxbase);
            z = z *(seqwordcount(proteinSample,'P')/maxbase)^(1/log(seqwordcount(proteinSample,'P')/maxbase));
            Z=[Z z]; 
        case('S')
            z = (seqwordcount(proteinSample,'S')/maxbase) * log(seqwordcount(proteinSample,'S')/maxbase);
            z = z *(seqwordcount(proteinSample,'S')/maxbase)^(1/log(seqwordcount(proteinSample,'S')/maxbase));
            Z=[Z z]; 
        case('T')
            z = (seqwordcount(proteinSample,'T')/maxbase) * log(seqwordcount(proteinSample,'T')/maxbase);
            z = z *(seqwordcount(proteinSample,'T')/maxbase)^(1/log(seqwordcount(proteinSample,'T')/maxbase));
            Z=[Z z]; 
        case('W')
            z = (seqwordcount(proteinSample,'W')/maxbase) * log(seqwordcount(proteinSample,'W')/maxbase);
            z = z *(seqwordcount(proteinSample,'W')/maxbase)^(1/log(seqwordcount(proteinSample,'W')/maxbase));
            Z=[Z z]; 
        case('Y')
            z = (seqwordcount(proteinSample,'Y')/maxbase) * log(seqwordcount(proteinSample,'Y')/maxbase);
            z = z *(seqwordcount(proteinSample,'Y')/maxbase)^(1/log(seqwordcount(proteinSample,'Y')/maxbase));
            Z=[Z z]; 
        case('V')
            z = (seqwordcount(proteinSample,'V')/maxbase) * log(seqwordcount(proteinSample,'V')/maxbase);
            z = z *(seqwordcount(proteinSample,'V')/maxbase)^(1/log(seqwordcount(proteinSample,'V')/maxbase));
            Z=[Z z];
        otherwise;
            break; 
    end
    Z = abs(Z)
end

%Calculation of locations and numbers of amino acids
indexA = find(proteinSample=='A')
indexR = find(proteinSample=='R')
indexN = find(proteinSample=='N')
indexD = find(proteinSample=='D')
indexC = find(proteinSample=='C')
indexQ = find(proteinSample=='Q')
indexE = find(proteinSample=='E')
indexG = find(proteinSample=='G')
indexH = find(proteinSample=='H')
indexI = find(proteinSample=='I')
indexL = find(proteinSample=='L')
indexK = find(proteinSample=='K')
indexM = find(proteinSample=='M')
indexF = find(proteinSample=='F')
indexP = find(proteinSample=='P')
indexS = find(proteinSample=='S')
indexT = find(proteinSample=='T')
indexW = find(proteinSample=='W')
indexY = find(proteinSample=='Y')
indexV = find(proteinSample=='V')
indexB = find(proteinSample=='B')
indexZ = find(proteinSample=='Z')
indexU = find(proteinSample=='U')



