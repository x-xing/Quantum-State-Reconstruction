function val=minor(matrix, rows, cols)
% Usage: m = minor(matrix, rows, cols)
%
% Calculates a minor of a given matrix, obtained by dropping all
% columns listed in "cols" and all rows listed in "rows" from
% matrix, then taking the determenent.  Matrix should be square and
% length(rows)==length(cols).

  matrix(rows,:) = [];
  matrix(:,cols) = [];
  val = det(matrix);
