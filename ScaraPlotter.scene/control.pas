// Constants
const
  NumJoints = 3;
  // Links
  l_1  = 0.6;   l_2  = 0.4;                 // length
  m_1  = 0.2;   m_2  = 0.2;   m_3  = 0.1;   // mass
  lc_1 = l_1/2; lc_2 = l_2/2;               // length to center of mass
  lx_1 = l_1;   lx_2 = l_2;   lx_3 = 0.02;  // dimensions
  ly_1 = 0.04;  ly_2 = 0.04;  ly_3 = 0.02;  //   (relative to frame oci_xci_yci)
  lz_1 = 0.04;  lz_2 = 0.04;  lz_3 = 0.05;
  // Motors
  n_1  = 1;     n_2  = 1;     n_3  = 1;     // gear reduction ratio
  km_1 = 0.56;  km_2 = 0.56;  km_3 = 0.56;  // torque constant
  jm_1 = 0;     jm_2 = 0;     jm_3 = 0;     // inertia
  bm_1 = 0.10;  bm_2 = 0.10;  bm_3 = 0.10;  // viscous constant

// Global Variables
var iZAxis, iJ1Axis, iJ2Axis, iPen: integer;



// Utilitary functions
// - set reference of  the joints (J1,J2: revolute; J3: prismatic)
procedure SetValuesJoints(Values: matrix);
begin
  SetAxisPosRef(0, iJ1Axis, Mgetv(Values, 0, 0));
  SetAxisPosRef(0, iJ2Axis, Mgetv(Values, 1, 0));
  SetAxisPosRef(0, iZAxis , Mgetv(Values, 2, 0));
end;

// - compute Denavit-Hartenberg (DH) matrix
function DHMat(a, alpha, d, theta: double): Matrix;
var ct, st, ca, sa: double;
    R: Matrix;
begin
  ct := cos(theta);
  st := sin(theta);
  ca := cos(alpha);
  sa := sin(alpha);

  R := Meye(4);
  MSetV(R,0,0,ct); MSetV(R,0,1,-st*ca); MSetV(R,0,2, st*sa); MSetV(R,0,3,a*ct);
  MSetV(R,1,0,st); MSetV(R,1,1, ct*ca); MSetV(R,1,2,-ct*sa); MSetV(R,1,3,a*st);
  MSetV(R,2,0, 0); MSetV(R,2,1, sa   ); MSetV(R,2,2, ca   ); MSetV(R,2,3,d   );
  result := R;
end;

// TODO: compute matrix D - inertia matrix
// TODO: compute diagonal matrix r^2 * Jm
// TODO: compute matrix M
// TODO: compute matrix C
// TODO: compute matrix B
// TODO: compute matrix Phi



// Forward Kinematics
function DK3(JointValues: matrix): matrix;
var
  A1, A2, A3: Matrix;
  P: Matrix;
  theta1, theta2, disp3: double;
begin
  theta1 := Mgetv(JointValues, 0, 0);
  theta2 := Mgetv(JointValues, 1, 0);
  disp3 := Mgetv(JointValues, 2, 0);

  A1 := DHMat(l_1,0,0,theta1);
  A2 := DHMat(l_2, rad(180), 0,theta2);
  A3 := DHMat(0,0,disp3,0);

  P := MMult(A1, A2);
  P := MMult(P, A3);

  result := P;
end;

// Inverse Kinematics
function IK3(XYZ: matrix): matrix;
var
  xc, yc, zc, c2: double;
  theta1, theta2, disp3: double;
begin
  xc := Mgetv(XYZ, 0, 0);
  yc := Mgetv(XYZ, 1, 0);
  zc := Mgetv(XYZ, 2, 0);

  c2 := (Power(xc,2) + power(yc,2) - power(l_1,2) - power(l_2,2))/(2*l_1*l_2);

  theta2 := ATan2(-sqrt(1-c2),c2);
  theta1 := ATan2(yc,xc) - ATan2(l_2*sin(theta2), l_1 + l_2*cos(theta2));
  disp3 := -zc;

  result := Mzeros(3, 1);
  MSetV(result, 0, 0, theta1);
  MSetV(result, 1, 0, theta2);
  MSetV(result, 2, 0, disp3);
end;



// Main Control Cycle: this procedure is called periodicaly 
// (default: 40 ms; can change it in Config > Control > Global > ScriptPeriod)
procedure Control;
var
  t: tcanvas;
  HTrans, XYZ, ORIENTATION, AuxJoints, ReqJoints, ReqValuesDK: matrix;
  U1, U2: double;
  PosPen: TPoint3D;

begin
  // Set color of the pen (RGB)
  if RCButtonPressed(1, 2) then
    SetSensorColor(0, 0, round(GetRCValue(1, 3)), round(GetRCValue(1, 4)), 
        round(GetRCValue(1, 5)));

  // Set position of the pen
  // - up
  if RCButtonPressed(3, 2) then
    SetAxisPosRef(0, iZAxis, GetRCValue(3, 3));
  // - down
  if RCButtonPressed(4, 2) then
    SetAxisPosRef(0, iZAxis, GetRCValue(4, 3));

  // Set angle of joints J1,J2
  if RCButtonPressed(6, 2) then begin
    SetAxisPosRef(0, iJ1Axis, Deg(GetRCValue(6, 3)));
    SetAxisPosRef(0, iJ2Axis, Deg(GetRCValue(6, 4)));
  end;
  if RCButtonPressed(7, 2) then begin
    SetAxisPosRef(0, iJ1Axis, 0);
    SetAxisPosRef(0, iJ2Axis, 0);
  end;

  // Set voltage of motor
  // (default controller should be disabled:
  //    Scene > Tag articulations > Tag controller > active='0' > Rebuild Scene)
  // - joint 1
  if RCButtonPressed(11, 2) then begin
    U1 := GetRCValue(11, 3);
    SetAxisVoltageRef(0, iJ1Axis, U1);
  end;
  // - joint 2
  if RCButtonPressed(12, 2) then begin
    U2 := GetRCValue(12, 3);
     SetAxisVoltageRef(0, iJ2Axis, U2);
  end;
  // - reset voltages
  if RCButtonPressed(11, 4) then begin
    U1 := 0;
    U2 := 0;
    SetAxisVoltageRef(0, iJ1Axis, U1);
    SetAxisVoltageRef(0, iJ1Axis, U2);
  end;



   // Forward Kinematics
  if RCButtonPressed(1, 9) then begin
    ReqValuesDK := Mzeros(3, 1);
    Msetv(ReqValuesDK, 0, 0, rad(GetRCValue(2, 9)));
    Msetv(ReqValuesDK, 1, 0, rad(GetRCValue(3, 9)));
    Msetv(ReqValuesDK, 2, 0, GetRCValue(4, 9));

    SetValuesJoints(ReqValuesDK);
    HTrans := DK3(ReqValuesDK);
    XYZ := Mzeros(3,1);

    MSetV(XYZ, 0, 0, Mgetv(HTrans, 0, 3));
    MSetV(XYZ, 1, 0, Mgetv(HTrans, 1, 3));
    MSetV(XYZ, 2, 0, Mgetv(HTrans, 2, 3));

    ORIENTATION := Meye(3);

    MSetV(ORIENTATION,0,0,Mgetv(HTrans,0,0)); MSetV(ORIENTATION,0,1,Mgetv(HTrans,0,1)); MSetV(ORIENTATION,0,2,Mgetv(HTrans,0,2));
    MSetV(ORIENTATION,1,0,Mgetv(HTrans,1,0)); MSetV(ORIENTATION,1,1,Mgetv(HTrans,1,1)); MSetV(ORIENTATION,1,2,Mgetv(HTrans,1,2));
    MSetV(ORIENTATION,2,0,Mgetv(HTrans,2,0)); MSetV(ORIENTATION,2,1,Mgetv(HTrans,2,1)); MSetV(ORIENTATION,2,2,Mgetv(HTrans,2,2));

    MatrixToRangeF(2, 10, XYZ, '%.3f');
    MatrixToRangeF(6, 7, ORIENTATION, '%.3f');
  end;

  // Inverse Kinematics
  if RCButtonPressed(1,14) then begin
    XYZ := Mzeros(3,1);

    Msetv(XYZ, 0, 0, GetRCValue(2, 14));
    Msetv(XYZ, 1, 0, GetRCValue(3, 14));
    Msetv(XYZ, 2, 0, GetRCValue(4, 14));

    ReqJoints := IK3(XYZ);
    SetValuesJoints(ReqJoints);
    AuxJoints := Mzeros(3,1);

    MSetV(AuxJoints, 0, 0, Deg(Mgetv(ReqJoints, 0, 0)));
    MSetV(AuxJoints, 1, 0, Deg(Mgetv(ReqJoints, 1, 0)));
    MSetV(AuxJoints, 2, 0, Mgetv(ReqJoints, 2, 0));

    MatrixToRangeF(2, 15, AuxJoints, '%.3f');
  end;

  SetRCValue(2, 8, format('%.3g',[Deg(GetAxisPos(0, iJ1Axis))]));
  SetRCValue(3, 8, format('%.3g',[Deg(GetAxisPos(0, iJ2Axis))]));
  SetRCValue(4, 8, format('%.3g',[GetAxisPos(0, iZAxis)]));

  PosPen := GetSolidPos(0, iPen);
  SetRCValue(2, 13, format('%.3g',[PosPen.x]));

  if GetRCValue(14, 2) = 0 then begin
    SetRCValue(3, 13, format('%.3g',[PosPen.y]));
    SetRCValue(4, 13, format('%.3g',[PosPen.z - 0.25]));
  end else begin
    SetRCValue(3, 13, format('%.3g',[PosPen.z]));
    SetRCValue(4, 13, format('%.3g',[-PosPen.y - 0.25]));
  end;

  t := GetSolidCanvas(0,0);
  t.brush.color := clwhite;
  t.textout(10,10, GetRCText(9, 2));
end;

// Initialization procedure: is called once when the script is started
procedure Initialize;
begin
  iJ1Axis := GetAxisIndex(0, 'Joint1', 0);
  iJ2Axis := GetAxisIndex(0, 'Joint2', 0);
  iZAxis  := GetAxisIndex(0, 'SlideZ', 0);
  iPen := GetSolidIndex(0, 'Pen');
end;
