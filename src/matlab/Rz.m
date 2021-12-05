function R = Rz(thz)
%RX

  R = [
    cos(thz) , -sin(thz) , 0 ;
    sin(thz) ,  cos(thz) , 0 ;
           0 ,         0 , 1 
  ];
end

