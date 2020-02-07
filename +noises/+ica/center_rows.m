function [CentMatrix, Mean] = center_rows(Matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mean = mean(Matrix, 2);
CentMatrix = Matrix - mean(Matrix, 2);
end