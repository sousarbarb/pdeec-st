function R = Rx(thx)
%RX

  R = [
    1 ,        0 ,         0 ;
    0 , cos(thx) , -sin(thx) ;
    0 , sin(thx) ,  cos(thx)
  ];
end

