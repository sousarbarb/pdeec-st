function [Aout] = Transpose(Ain)
%TRANSPOSE

  [m,n] = size(Ain);
  %Aout = Ain;
  Aout = vpa(zeros(n,m));
  for i=1:m
    for j=1:n
      Aout(j,i) = Ain(i,j);
    end
  end
end

