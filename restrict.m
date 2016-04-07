function [v_new,row_interweave] = restrict( v_old)
%RESTRICT - performs a restriction operation from the fine grid to the
%coarse grid. Type can be injection or full weighting. 
%Types: 
%           fw - full weighting
%           in - injection


        %FULL WEIGHTING - does full weighting on the fine grid
        %   Set up restriction matrix
        %   Set up matrix of (n/2-1)x(n-1)
        %   Form (n/2-1)x(n/2-1) and n/2x(n/2-1)
        %   Take the TRANSPOSE of the INTERPOLATION matrix 
        
        n = length(v_old); %v_old is length n-1
        n = n/2 -1/2;       %adjust length so it is n/2-1
        d = 2*ones(n,1);
        d1 = ones(n+1,1);
        
        %create matrix with twos on diagonal size n/2-1 by n/2-1
        two_matrix = spdiags(d, 0, n+1,n);
        %create matrix with twos on diagonal size n/2 by n/2-1
        ones_matrix = spdiags([d1 d1],-1:0, n+1,n);
       % ones_matrix = full(ones_matrix);
        rows = size(two_matrix,1) + size(ones_matrix,1); %adding the number of rows of the two matrices
        %two_matrix = full(two_matrix);
        row_interweave = reshape([ones_matrix(:) two_matrix(:) ]',rows, []);
        row_restrict = row_interweave(1:end-1,:)';
        
        %row_restrict
        %if v_old is a row vector, put your thang down flip it and reverse it
        if iscolumn(v_old) == 0
            v_old = v_old';
        end
        
        v_new = 1/4*row_restrict*v_old;
        %v_new has length n-1


end

