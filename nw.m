function alignedSeq = nw(x)
    Dataset1 = char(x(1));
    Dataset2 = char(x(2));
    
% Calculate length of sequences 
% gap penalities is %1
Seq1len = length(Dataset1)+1;
Seq2len = length(Dataset2)+1;

% Initilization
% Trackback 
% Matrix created
NeedlemanWunsch = zeros(Seq1len, Seq2len);
NeedlemanWunschSolution = zeros(Seq1len-1, Seq2len-1);

    %initialize scores
    a=1;
    b=2;
    c=3;


    for n = 2:Seq1len
        NeedlemanWunsch(n,1) = (-n)*b;
    end

    for n = 2:Seq2len
        NeedlemanWunsch(1,n) = (-n)*b;
    end

    for i = 2:Seq1len
        d1 = Dataset1(i-1);
        for j = 2:Seq2len
            d2 = Dataset2(j-1);

             if d1 == d2
                score = a;
            else 
                score = -c;
            end

            % Optimal score
            scoreCase = zeros(1,3);
            scoreCase(1) = NeedlemanWunsch(i-1,j-1) + score;
            scoreCase(2) = NeedlemanWunsch(i-1,j) - b;
            scoreCase(3) = NeedlemanWunsch(i,j-1) - b;

            % Matrix now shows the score and direction
           [score,Direction] = max(scoreCase);
            NeedlemanWunsch(i,j)=score;
            NeedlemanWunschSolution(i-1,j-1)=Direction;
        end
    end

    % Sequences are built
    alignedSeq1 = '';
    alignedSeq2 = '';
    i = Seq1len-1;
    j = Seq2len-1;
    
    % Traceback
   while j > 0 || i > 0
        if i>0 && j>0
            
            % Diagonal
           if NeedlemanWunschSolution(i,j)==1   
               alignedSeq1 = strcat(Dataset1(i),alignedSeq1);
               alignedSeq2 = strcat(Dataset2(j),alignedSeq2);
               i=i-1;
               j=j-1;
               continue
           end
           
           % Up
           if NeedlemanWunschSolution(i,j)==2   
               alignedSeq1 = strcat(Dataset1(i),alignedSeq1);
               alignedSeq2 = strcat('-',alignedSeq2);
               i=i-1;
               continue
           end
           
           % left
           if NeedlemanWunschSolution(i,j)==3   
               alignedSeq1 = strcat('-',alignedSeq1);
               alignedSeq2 = strcat(Dataset2(j),alignedSeq2);
               j=j-1;
               continue
           end
        else
            break;
        end
    end
    
    alignedSeq = {alignedSeq1 , alignedSeq2};
end
