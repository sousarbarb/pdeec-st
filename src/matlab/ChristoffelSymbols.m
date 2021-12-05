function [cijk] = ChristoffelSymbols(D,q,i,j,k)
%CHRISTOFFELSYMBOLS

  cijk = ( diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)) ) / 2;
end

