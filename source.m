% Read the files
Input1 = fastaread('Human_HOX.fa');
Input2 = fastaread('Fly_HOX.fa');
MatchFile = fopen('Match_HOX.txt', 'r');

% Create arrays after the match file is read 
if MatchFile>-1
    MatchTable = fscanf(MatchFile, '%d');
    i=1;
    s=1;
    finalSeq1 = '';
    finalSeq2 = '';
    
    while i<length(MatchTable)

        Input1Matches(s) = MatchTable(i);
        Input1Matches(s+1) = MatchTable(i+1);
        Input2Matches(s) = MatchTable(i+2);
        Input2Matches(s+1) = MatchTable(i+3);
        s=s+2;
        i=i+4;
    end
    
   i=0;
    
   while i<length(Input1Matches)
        
        % Obtain the sequence index from the match file
       if i == 0
            x1 = 1;
            y1 = 1;
        else
            x1 = Input1Matches(i);
            y1 = Input2Matches(i);
        end
        x2 = Input1Matches(i+1)-1;
        y2 = Input2Matches(i+1)-1;

        Seq1 = Input1.Sequence(x1:x2);
        Seq2 = Input2.Sequence(y1:y2);

        % Call to the NeedlemanWunsch algorithm 
        BothSeq = {Seq1,Seq2};
        alignedSeq = nw(BothSeq);
        
        % The sequences are put together
        notedMatch1 = Input1.Sequence(Input1Matches(i+1):Input1Matches(i+2)-1);
        notedMatch2 = Input2.Sequence(Input2Matches(i+1):Input2Matches(i+2)-1);
        finalSeq1 = strcat(finalSeq1,char(alignedSeq(1)));  
        finalSeq1 = strcat(finalSeq1,notedMatch1);          
        finalSeq2 = strcat(finalSeq2,char(alignedSeq(2)));  
        finalSeq2 = strcat(finalSeq2,notedMatch2);          
        
        fprintf('Alignment = %s\n',char(alignedSeq(1)));
        
        % Score computation 
         score = 0;
        alignedSeq1 = char(alignedSeq(1));
        alignedSeq2 = char(alignedSeq(2));
        for (s = 1:length(alignedSeq1))
            if (alignedSeq1(s)=='-' || alignedSeq2(s)=='-')
                score = score - 2;
            elseif (alignedSeq1(s)==alignedSeq2(s))
                score = score + 1;
            else
                score = score - 3;
            end                
        end
        fprintf('Score = %i\n\n',score);
            
        
        i=i+2;
    end
    
    % The end of the sequences are aligned.
    Seq1 = Input1.Sequence(Input1Matches(i):length(Input1.Sequence));
    Seq2 = Input2.Sequence(Input2Matches(i):length(Input2.Sequence));

    % Call to the NeedlemanWunsch algorithm
    BBothSeq = {Seq1,Seq2};
    alignedSeq = nw(BothSeq);
    finalSeq1 = strcat(finalSeq1,char(alignedSeq(1))); 
    finalSeq2 = strcat(finalSeq2,char(alignedSeq(2)));  
    
    % Score computation 
        score = 0;
        alignedSeq1 = char(alignedSeq(1));
        alignedSeq2 = char(alignedSeq(2));
        for (s = 1:length(alignedSeq1))
            if (alignedSeq1(s)=='-' || alignedSeq2(s)=='-')
                score = score - 2;
            elseif (alignedSeq1(s)==alignedSeq2(s))
                score = score + 1;
            else
                score = score - 3;
            end                
        end
        
        fprintf('Alignment = %s\n',char(alignedSeq(1)));
        fprintf('Aignment = %s\n',char(alignedSeq(2)));
        fprintf('Score = %i\n\n',score);
    
    fprintf('Anchor Alignment = %s\n',char(finalSeq1));
    fprintf('Anchor Alignment = %s\n',char(finalSeq2));
    
    
else
    % Entire sequence is used for computation
    BothSeq = {Input1.Sequence,Input2.Sequence};
    alignedSeq = nw(BothSeq);
    
    fprintf('Global Alignment = %s\n',char(alignedSeq(1)));
    fprintf('Global Alignment = %s\n',char(alignedSeq(2)));
    
    % Score computation
    score = 0;
    alignedSeq1 = char(alignedSeq(1));
    alignedSeq2 = char(alignedSeq(2));
    for (s = 1:length(alignedSeq1))
        if (alignedSeq1(s)=='-' || alignedSeq2(s)=='-')
            score = score - 2;
        elseif (alignedSeq1(s)==alignedSeq2(s))
            score = score + 1;
        else
            score = score - 3;
        end                
    end
    fprintf('Score = %i\n\n',score);
    
end

%% part 4: Repeating alignment 10,000 times

Input3 = fastaread('Human_PAX.fa');
Input4 = fastaread('Fly_PAX.fa');

permSeq1 = Input3.Sequence;
permSeq2 = Input4.Sequence;
r=1;
index1 = randperm(length(permSeq1));
index2 = randperm(length(permSeq2));
for i = 1:10000
   
    % Pick a sequence to permute 
    i
    if rand(1)==0
        
        % Random selection of an index from the permuted ones
        r = round(rand(1)*(length(permSeq1)-1));
        r = r + 1;

        % Locate two postions and exchange them
       x = double(index1(r));
        if (r==length(permSeq1))    r = 1;
        end
        x2 = double(index1(r+1));
        temp = permSeq1(x);
        permSeq1(x)= permSeq1(x2);
        permSeq1(x2) = temp;  
    else       
        % Random selection of an index from the permuted ones
        r = round(rand(1)*(length(permSeq2)-1));
        r = r + 1;
        
        % Locate two postions and exchange them
        x = double(index2(r));
        if (r==length(permSeq2))    r = 1;
        end
        x2 = double(index2(r+1));
        temp = permSeq2(x);
        permSeq2(x)= permSeq2(x2);
        permSeq2(x2) = temp;
    end
    
    % Using NeedlemanWunsch, compute the sequences 
    BothSeq = {permSeq1,permSeq2};
    alignedSeq = nw(BothSeq);
    
    fprintf('Global Alignment = %s\n',char(alignedSeq(1)));
    fprintf('Global Alignment = %s\n',char(alignedSeq(2)));
    
    % Score computation
    score = 0;
    alignedSeq1 = char(alignedSeq(1));
    alignedSeq2 = char(alignedSeq(2));
    for (s = 1:length(alignedSeq1))
        if (alignedSeq1(s)=='-' || alignedSeq2(s)=='-')
            score = score - 2;
        elseif (alignedSeq1(s)==alignedSeq2(s))
            score = score + 1;
        else
            score = score - 3;
        end                
    end
    totalHOXScore(i) = score;
    fprintf('Score = %i\n\n',score);
end


permSeq1 = Input1.Sequence;
permSeq2 = Input2.Sequence;
r=1;
index1 = randperm(length(permSeq1));
index2 = randperm(length(permSeq2));
for i = 1:10000
   
    % Pick a sequence to permute 
    i
    if rand(1)==0
        
        % Random selection of an index from the permuted ones
        r = round(rand(1)*(length(permSeq1)-1));
        r = r + 1;
        
        
        % Locate two postions and exchange them
        if (r==length(permSeq1))    r = 1;
        end
        x2 = double(index1(r+1));
        temp = permSeq1(x);
        permSeq1(x)= permSeq1(x2);
        permSeq1(x2) = temp;      
    else     
         % Random selection of an index from the permuted ones
        r = round(rand(1)*(length(permSeq2)-1));
        r = r + 1;
        
        % Locate two postions and exchange them
        x = double(index2(r));
        if (r==length(permSeq2))    r = 1;
        end
        x2 = double(index2(r+1));
        temp = permSeq2(x);
        permSeq2(x)= permSeq2(x2);
        permSeq2(x2) = temp;
    end
    
    % Using NeedlemanWunsch, compute the sequences 
    BothSeq = {permSeq1,permSeq2};
    alignedSeq = nw(BothSeq);
    
    fprintf('Global Alignment = %s\n',char(alignedSeq(1)));
    fprintf('Global Alignment = %s\n',char(alignedSeq(2)));
    
    % Score computation
    score = 0;
    alignedSeq1 = char(alignedSeq(1));
    alignedSeq2 = char(alignedSeq(2));
    for (s = 1:length(alignedSeq1))
        if (alignedSeq1(s)=='-' || alignedSeq2(s)=='-')
            score = score - 2;
        elseif (alignedSeq1(s)==alignedSeq2(s))
            score = score + 1;
        else
            score = score - 3;
        end                
    end
    totalPAXScore(i) = score;
    fprintf('Score = %i\n\n',score);
end