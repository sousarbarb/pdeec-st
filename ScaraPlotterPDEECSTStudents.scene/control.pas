{******************************************************************************
 * LABORATORIAL WORK X - INDUSTRIAL ROBOTICS - MIEEC AT FEUP
 * ANTÓNIO PAULO MOREIRA, PAULO G. COSTA, RICARDO B. SOUSA
 ******************************************************************************}

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
  rho_1 = m_1/(lx_1*ly_1*lz_1);             // mass density
  rho_2 = m_2/(lx_2*ly_2*lz_2);
  rho_3 = m_3/(lx_3*ly_3*lz_3);
  // Motors
  n_1  = 1;     n_2  = 1;     n_3  = 1;     // gear reduction ratio
  km_1 = 0.56;  km_2 = 0.56;  km_3 = 0.56;  // torque constant
  jm_1 = 0;     jm_2 = 0;     jm_3 = 0;     // inertia
  bm_1 = 0.10;  bm_2 = 0.10;  bm_3 = 0.10;  // viscous constant
  ri_1 = 0.12;  ri_2 = 0.12;  ri_3 = 0.12;  // internal resistance
  // Gravity
  gacc = 9.8;                               // gravity acceleration (m/s^2)

// Global Variables
var iZAxis, iJ1Axis, iJ2Axis, iPen: integer;
    Q_meas, Qd1_meas, W_mat, K0_mat, K1_mat: Matrix;
    Q_desir, Qd1_desir, Qd2_desir, XYZ_desir: Matrix;
    InvDynON: boolean;



{******************************************************************************
 * UTILITARY FUNCTIONS
 * - SetValuesJoints(Values: matrix):
 *     set the reference of a pre-defined and activated PID controller for
 *     each joint
 *
 *     [in] Values 3x1 matrix containing the desired reference
 *
 * - DHMat(theta, d, a, alpha: double): Matrix
 *     compute Denavit-Hartenberg (DH) matrix
 *
 *     [in] theta,d,a,alpha parameters defined in the DH convention
 *
 * - RotZ(theta: double): Matrix
 *     compute a rotation matrix in z-axis
 *
 *     [in] theta angle (rad) of the rotation
 *
 * - SetJointsControllerState(activate: boolean)
 *     change the controllers state. Only with the controllers deactivated it
 *     is possible to control the joints on voltage or torque
 *
 *     [in] activate new state desired for the controllers
 *
 * - SetTorque(Torque: matrix)
 *     set the torque value for each joint (the controllers must be
 *     deactivated)
 *
 *     [in] Torque 3x1 matrix containing the desired torque for each joint
 ******************************************************************************}
procedure SetValuesJoints(Values: matrix);
begin
  // Controller mode: position control
  SetMotorControllerMode(0, iJ1Axis, 'pidposition');
  SetMotorControllerMode(0, iJ2Axis, 'pidposition');
  SetMotorControllerMode(0, iZAxis , 'pidposition');
  // Controller parameters:          ki  kp  kp  kf
  SetMotorControllerPars(0, iJ1Axis,  0, 25,  3,  0);
  SetMotorControllerPars(0, iJ2Axis,  0, 50,  4,  0);
  SetMotorControllerPars(0, iZAxis ,  0,  1,  5,  0);
  // Set position reference
  SetAxisPosRef(0, iJ1Axis, Mgetv(Values, 0, 0));
  SetAxisPosRef(0, iJ2Axis, Mgetv(Values, 1, 0));
  SetAxisPosRef(0, iZAxis , Mgetv(Values, 2, 0));
end;

function DHMat(theta, d, a, alpha: double): Matrix;
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

function RotZ(theta: double): Matrix;
var ct, st: double;
    R: Matrix;
begin
  ct := cos(theta);
  st := sin(theta);
  R := Mzeros(3,3);

  MSetV(R,0,0,ct); MSetV(R,0,1,-st);
  MSetV(R,1,0,st); MSetV(R,1,1, ct);
                                     MSetV(R,2,2,1);

  result := R;
end;

procedure SetJointsControllerState(activate: boolean);
begin
  SetMotorControllerState(0, iJ1Axis, activate);
  SetMotorControllerState(0, iJ2Axis, activate);
  SetMotorControllerState(0, iZAxis , activate);
end;

procedure SetTorque(Torque: matrix);
begin
  SetAxisTorqueRef(0, iJ1Axis, Mgetv(Torque, 0, 0));
  SetAxisTorqueRef(0, iJ2Axis, Mgetv(Torque, 1, 0));
  SetAxisTorqueRef(0, iZAxis , Mgetv(Torque, 2, 0));
end;



{******************************************************************************
 * INVERSE DYNAMICS (ID)
 * - SetIDParameters
 *     set the controller gains K0 and K1 required for the outer loop
 *
 * - IDDMat(_Q: Matrix): Matrix
 *     compute the inertia matrix D of the SCARA manipulator
 *
 *     [in] _Q 3x1 matrix measured position of each joint
 *
 * - IDJmMat: Matrix
 *     compute the diagonal matrix r^2 * Jm
 *
 * - IDMMat(_Q: Matrix): Matrix
 *     compute the matrix M = D + r^2 * Jm
 *
 *     [in] _Q 3x1 matrix measured position of each joint
 *
 * - IDCMat(_Q,_Qd1: Matrix): Matrix
 *     compute the matrix C
 *
 *     [in] _Q   3x1 matrix measured position of each joint
 *     [in] _Qd1 3x1 matrix measured velocity of each joint
 *
 * - IDBMat: Matrix
 *     compute the diagonal matrix Bm + KbKm/R
 *
 * - IDPhiMatHorizontalSCARA(_Q: Matrix): Matrix
 *     compute the gravity compensation for the Horizontal SCARA manipulator
 *
 *     [in] _Q 3x1 matrix measured position of each joint
 *
 * - IDPhiMatVerticalSCARA(_Q: Matrix): Matrix
 *     compute the gravity compensation for the Vertical SCARA manipulator
 *
 *     [in] _Q 3x1 matrix measured position of each joint
 ******************************************************************************}
procedure SetIDParameters;
begin
  // wn
  W_mat := RangeToMatrix(12,10,3,1);
  // K0 = wn ^ 2
  K0_mat := Mzeros(3,3);
  MSetV(K0_mat,0,0,MGetV(W_mat,0,0)*MGetV(W_mat,0,0));
  MSetV(K0_mat,1,1,MGetV(W_mat,1,0)*MGetV(W_mat,1,0));
  MSetV(K0_mat,2,2,MGetV(W_mat,2,0)*MGetV(W_mat,2,0));
  // K1 = 2 * wn
  K1_mat := Mzeros(3,3);
  MSetV(K1_mat,0,0,2*MGetV(W_mat,0,0));
  MSetV(K1_mat,1,1,2*MGetV(W_mat,1,0));
  MSetV(K1_mat,2,2,2*MGetV(W_mat,2,0));
end;

function IDDMat(_Q: Matrix): Matrix;
var q1_, q2_, q3_: double;
    cq1, cq2, cq1pq2, sq1, sq2, sq1pq2: double;
    Jvc_1, Jvc_2, Jvc_3: Matrix;
    Jwc_1, Jwc_2, Jwc_3: Matrix;
    Rc_1, Rc_2, Rc_3: Matrix;
    I_1, I_2, I_3: Matrix;
begin
  q1_ := MGetV(_Q, 0, 0);
  q2_ := MGetV(_Q, 1, 0);
  q3_ := MGetV(_Q, 2, 0);
  cq1 := cos(q1_);
  sq1 := sin(q1_);
  cq2 := cos(q2_);
  sq2 := sin(q2_);
  cq1pq2 := cos(q1_+q2_);
  sq1pq2 := sin(q1_+q2_);

  // Matrices Jvc_i (linear speed jacobian relative to oci_xci_yci_zci)
  Jvc_1 := Mzeros(3,3);
  Msetv(Jvc_1,0,0,-lc_1*sq1);
  Msetv(Jvc_1,1,0, lc_1*cq1);
  
  Jvc_2 := Mzeros(3,3);
  Msetv(Jvc_2,0,0,-l_1*sq1-lc_2*sq1pq2); Msetv(Jvc_2,0,1,-lc_2*sq1pq2);
  Msetv(Jvc_2,1,0, l_1*cq1+lc_2*cq1pq2); Msetv(Jvc_2,1,1, lc_2*cq1pq2);
  
  Jvc_3 := Mzeros(3,3);
  Msetv(Jvc_3,0,0,-l_1*sq1-l_2*sq1pq2); Msetv(Jvc_3,0,1,-l_2*sq1pq2);
  Msetv(Jvc_3,1,0, l_1*cq1+l_2*cq1pq2); Msetv(Jvc_3,1,1, l_2*cq1pq2);
                                                                      Msetv(Jvc_3,2,2,-1);

  // Matrices Jwc_i (angular speed jacobian relative to oci_xci_yci_zci)
  Jwc_1 := Mzeros(3,3);
  Msetv(Jwc_1,2,0,1);

  Jwc_2 := Mzeros(3,3);
  Msetv(Jwc_2,2,0,1); Msetv(Jwc_2,2,1,1);

  Jwc_3 := Jwc_2;

  // Matrices Rc_i (orientation matrices relative to oci_xci_yci_zci)
  Rc_1 := RotZ(q1_);
  
  Rc_2 := Mzeros(3,3);
  MSetV(Rc_2,0,0,cq1pq2); MSetV(Rc_2,0,1, sq1pq2);
  MSetV(Rc_2,1,0,sq1pq2); MSetV(Rc_2,1,1,-cq1pq2);
                                                   MSetV(Rc_2,2,2,-1);
  
  Rc_3 := Rc_2;

  // Matrices I_i (inertia tensors relative to oci_xci_yci_zci)
  I_1 := Mzeros(3,3);
  MSetV(I_1,0,0,rho_1*lx_1*ly_1*lz_1*(ly_1*ly_1+lz_1*lz_1)/12);
  MSetV(I_1,1,1,rho_1*lx_1*ly_1*lz_1*(lx_1*lx_1+lz_1*lz_1)/12);
  MSetV(I_1,2,2,rho_1*lx_1*ly_1*lz_1*(ly_1*ly_1+lx_1*lx_1)/12);

  I_2 := Mzeros(3,3);
  MSetV(I_2,0,0,rho_2*lx_2*ly_2*lz_2*(ly_2*ly_2+lz_2*lz_2)/12);
  MSetV(I_2,1,1,rho_2*lx_2*ly_2*lz_2*(lx_2*lx_2+lz_2*lz_2)/12);
  MSetV(I_2,2,2,rho_2*lx_2*ly_2*lz_2*(ly_2*ly_2+lx_2*lx_2)/12);

  I_3 := Mzeros(3,3);
  MSetV(I_3,0,0,rho_3*lx_3*ly_3*lz_3*(ly_3*ly_3+lz_3*lz_3)/12);
  MSetV(I_3,1,1,rho_3*lx_3*ly_3*lz_3*(lx_3*lx_3+lz_3*lz_3)/12);
  MSetV(I_3,2,2,rho_3*lx_3*ly_3*lz_3*(ly_3*ly_3+lx_3*lx_3)/12);

  // Inertia Matrix D(Q) (Q = [q1 q2 q3]')
  result := MMultReal( MMult(MTran(Jvc_1),Jvc_1) , m_1 );
  result := MAdd( result , MMultReal( MMult(MTran(Jvc_2),Jvc_2) , m_2 ) );
  result := MAdd( result , MMultReal( MMult(MTran(Jvc_3),Jvc_3) , m_3 ) );
  result := MAdd( result , MMult(MMult(MMult(MTran(Jwc_1),Rc_1),I_1),MMult(MTran(Rc_1),Jwc_1)) );
  result := MAdd( result , MMult(MMult(MMult(MTran(Jwc_2),Rc_2),I_2),MMult(MTran(Rc_2),Jwc_2)) );
  result := MAdd( result , MMult(MMult(MMult(MTran(Jwc_3),Rc_3),I_3),MMult(MTran(Rc_3),Jwc_3)) );
end;

function IDJmMat: Matrix;
begin
  result := Mzeros(3,3);
  MSetV(result,0,0,n_1*n_1*jm_1);
  MSetV(result,1,1,n_2*n_2*jm_2);
  MSetV(result,2,2,n_3*n_3*jm_3);
end;

function IDMMat(_Q: Matrix): Matrix;
begin
  result := MAdd( IDDMat(_Q) , IDJmMat );
end;

function IDCMat(_Q,_Qd1: Matrix): Matrix;
var q1_, q2_, q3_: double;
    q1d1_, q2d1_, q3d1_: double;
    cq1, cq2, cq1pq2, sq1, sq2, sq1pq2: double;
    c211, c112, c121, c221: double;
begin
  q1_ := MGetV(_Q, 0, 0);
  q2_ := MGetV(_Q, 1, 0);
  q3_ := MGetV(_Q, 2, 0);
  q1d1_ := MGetV(_Qd1, 0, 0);
  q2d1_ := MGetV(_Qd1, 1, 0);
  q3d1_ := MGetV(_Qd1, 2, 0);
  cq1 := cos(q1_);
  sq1 := sin(q1_);
  cq2 := cos(q2_);
  sq2 := sin(q2_);
  cq1pq2 := cos(q1_+q2_);
  sq1pq2 := sin(q1_+q2_);

  // Christoffel symbols
  c211 := m_3 * cq1pq2 * l_2  * ( l_1 * sq1 + l_2  * sq1pq2 ) -
          m_2 * sq1pq2 * lc_2 * ( l_1 * cq1 + lc_2 * cq1pq2 ) -
          m_3 * sq1pq2 * l_2  * ( l_1 * cq1 + l_2  * cq1pq2 ) +
          m_2 * cq1pq2 * lc_2 * ( l_1 * sq1 + lc_2 * sq1pq2 );
  c121 := c211;
  c221 := c211;
  c112 := m_3 * sq1pq2 * l_2  * ( l_1 * cq1 + l_2  * cq1pq2 ) +
          m_2 * sq1pq2 * lc_2 * ( l_1 * cq1 + lc_2 * cq1pq2 ) -
          m_3 * cq1pq2 * l_2  * ( l_1 * sq1 + l_2  * sq1pq2 ) -
          m_2 * cq1pq2 * lc_2 * ( l_1 * sq1 + lc_2 * sq1pq2 );

  // Matrix C
  result := Mzeros(3,3);
  MSetV(result,0,0,c211*q2d1_); MSetV(result,0,1,c121*q1d1_+c221*q2d1_);
  MSetV(result,1,0,c112*q1d1_);
end;

function IDBMat: Matrix;
begin
  result := Mzeros(3,3);
  MSetV(result,0,0,bm_1+km_1*km_1/ri_1);
  MSetV(result,1,1,bm_2+km_2*km_2/ri_2);
  MSetV(result,2,2,bm_3+km_3*km_3/ri_3);
end;

function IDPhiMatHorizontalSCARA(_Q: Matrix): Matrix;
begin
  result := Mzeros(3,1);
  MSetV(result,2,0,-gacc*m_3);
end;

function IDPhiMatVerticalSCARA(_Q: Matrix): Matrix;
var q1_, q2_, q3_: double;
    cq1, cq2, cq1pq2, sq1, sq2, sq1pq2: double;
begin
  q1_ := MGetV(_Q, 0, 0);
  q2_ := MGetV(_Q, 1, 0);
  q3_ := MGetV(_Q, 2, 0);
  cq1 := cos(q1_);
  cq1pq2 := cos(q1_+q2_);

  result := Mzeros(3,1);
  MSetV(result,0,0,gacc*lc_1*cq1*m_1 +
                   gacc*l_1 *cq1*m_2 + gacc*lc_2*cq1pq2*m_2+
                   gacc*l_1 *cq1*m_3 + gacc*l_2 *cq1pq2*m_3);
  MSetV(result,1,0,gacc*lc_2*cq1pq2*m_2 +
                   gacc*l_2 *cq1pq2*m_3);
  MSetV(result,2,0,0);
end;



{******************************************************************************
 * KINEMATICS
 * - DK3(_Q: Matrix): Matrix
 *     compute the forward kinematics of the SCARA manipulator
 *
 *     [in] _Q 3x1 matrix measured position of each joint
 *
 * - IK3(XYZ: Matrix): Matrix
 *     compute the inverse kinematics of the SCARA manipulator
 *
 *     [in] XYZ 3x1 matrix desired position of the end-effector
 ******************************************************************************}
function DK3(_Q: Matrix): Matrix;
var A1, A2, A3: Matrix;
    P: Matrix;
    q1_, q2_, q3_: double;
begin
  q1_ := MGetV(_Q, 0, 0);
  q2_ := MGetV(_Q, 1, 0);
  q3_ := MGetV(_Q, 2, 0);

  A1 := DHMat(q1_,0  ,l_1,0       );
  A2 := DHMat(q2_,0  ,l_2,rad(180));
  A3 := DHMat(0  ,q3_,0  ,0       );

  P := MMult(A1,A2);
  P := MMult(P,A3);

  result := P;
end;

function IK3(XYZ: Matrix): Matrix;
var
  xc, yc, zc, c2: double;
  q1_, q2_, q3_: double;
begin
  xc := Mgetv(XYZ, 0, 0);
  yc := Mgetv(XYZ, 1, 0);
  zc := Mgetv(XYZ, 2, 0);

  c2 := (Power(xc,2) + power(yc,2) - power(l_1,2) - power(l_2,2))/(2*l_1*l_2);

  q2_ := ATan2(-sqrt(1-power(c2,2)),c2);
  q1_ := ATan2(yc,xc) - ATan2(l_2*sin(q2_), l_1 + l_2*cos(q2_));
  q3_ := -zc;

  result := Mzeros(3, 1);
  MSetV(result, 0, 0, q1_);
  MSetV(result, 1, 0, q2_);
  MSetV(result, 2, 0, q3_);
end;



{******************************************************************************
 * CONTROL
 * Main Control Cycle: this procedure is called periodically
 * (default: 40 ms; can change it in Config > Control > Global > ScriptPeriod)
 ******************************************************************************}
procedure Control;
var
  t: tcanvas;
  HTrans, XYZ, ORIENTATION, AuxJoints, ReqJoints, ReqValuesDK: matrix;
  U1, U2, U3, T1, T2, T3: double;
  PosPen: TPoint3D;
  M_mat, C_mat, B_mat, Phi_mat: Matrix;
  Aq_mat, R_mat, U_mat: Matrix;
  Err_mat: Matrix;

begin
  // Initialization
  // - joint position
  MSetV(Q_meas,0,0, GetAxisPos(0, iJ1Axis) );
  MSetV(Q_meas,1,0, GetAxisPos(0, iJ2Axis) );
  MSetV(Q_meas,2,0, GetAxisPos(0, iZAxis ) );
  // - joint speed
  MSetV(Qd1_meas,0,0, GetAxisSpeed(0, iJ1Axis) );
  MSetV(Qd1_meas,1,0, GetAxisSpeed(0, iJ2Axis) );
  MSetV(Qd1_meas,2,0, GetAxisSpeed(0, iZAxis ) );
  // - pen position
  PosPen := GetSolidPos(0, iPen);
  // - canvas
  t := GetSolidCanvas(0,0);
  t.brush.color := clwhite;
  t.textout(10,10, GetRCText(8, 2));



  // Set color of the pen (RGB)
  if RCButtonPressed(1, 2) then
    SetSensorColor(0, 0, round(GetRCValue(1, 3)), round(GetRCValue(1, 4)), 
        round(GetRCValue(1, 5)));

  // Set position of the pen
  // - up
  if RCButtonPressed(3, 2) then begin
    SetJointsControllerState(true);
    SetAxisPosRef(0, iZAxis, GetRCValue(3, 3));
  end;
  // - down
  if RCButtonPressed(4, 2) then begin
    SetJointsControllerState(true);
    SetAxisPosRef(0, iZAxis, GetRCValue(4, 3));
  end;

  // Reset angle of joints J1,J2
  if RCButtonPressed(6, 2) then begin
    SetJointsControllerState(true);
    SetValuesJoints(MZeros(3,1));
  end;

  // Set voltage of motor
  // (default controller should be disabled:
  //    Scene > Tag articulations > Tag controller > active='0' > Rebuild Scene)
  // - joint 1
  if RCButtonPressed(12, 2) then begin
    U1 := GetRCValue(12, 3);
    SetJointsControllerState(false);
    SetAxisVoltageRef(0, iJ1Axis, U1);
  end;
  // - joint 2
  if RCButtonPressed(13, 2) then begin
    U2 := GetRCValue(13, 3);
    SetJointsControllerState(false);
    SetAxisVoltageRef(0, iJ2Axis, U2);
  end;
  // - joint 3
  if RCButtonPressed(14, 2) then begin
    U3 := GetRCValue(14, 3);
    SetJointsControllerState(false);
    SetAxisVoltageRef(0, iZAxis , U3);
  end;
  // - reset voltages
  if RCButtonPressed(12, 4) then begin
    U1 := 0;
    U2 := 0;
    U3 := 0;
    SetJointsControllerState(false);
    SetAxisVoltageRef(0, iJ1Axis, U1);
    SetAxisVoltageRef(0, iJ2Axis, U2);
    SetAxisVoltageRef(0, iZAxis , U3);
  end;

  // Set torque of motor
  // (default controller should be disabled:
  //    Scene > Tag articulations > Tag controller > active='0' > Rebuild Scene)
  // - joint 1
  if RCButtonPressed(16, 2) then begin
    T1 := GetRCValue(16, 3);
    SetJointsControllerState(false);
    SetAxisTorqueRef(0, iJ1Axis, T1);
  end;
  // - joint 2
  if RCButtonPressed(17, 2) then begin
    T2 := GetRCValue(17, 3);
    SetJointsControllerState(false);
    SetAxisTorqueRef(0, iJ2Axis, T2);
  end;
  // - joint 3
  if RCButtonPressed(18, 2) then begin
    T3 := GetRCValue(18, 3);
    SetJointsControllerState(false);
    SetAxisTorqueRef(0, iZAxis , T3);
  end;
  // - reset voltages
  if RCButtonPressed(16, 4) then begin
    T1 := 0;
    T2 := 0;
    T3 := 0;
    SetJointsControllerState(false);
    SetAxisTorqueRef(0, iJ1Axis, T1);
    SetAxisTorqueRef(0, iJ2Axis, T2);
    SetAxisTorqueRef(0, iZAxis , T3);
  end;



  // Forward Kinematics
  if RCButtonPressed(1, 9) then begin
    SetJointsControllerState(true);
    
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
    SetJointsControllerState(true);

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



  // Inverse Dynamics
  if RCButtonPressed(10,7) then begin
    if InvDynON then begin
      InvDynON := false;
      SetJointsControllerState(true);
    end else begin
      InvDynON := true;
      SetJointsControllerState(false);
      SetIDParameters;
      // Requested set point
      XYZ_desir := RangeToMatrix(12,9,3,1);
      Q_desir := IK3(XYZ_desir);
      Qd1_desir := Mzeros(3,1);
      Qd2_desir := Mzeros(3,1);
    end;
  end;
  if InvDynON AND RCButtonPressed(10,9) then begin
      SetIDParameters;
      // Requested set point
      XYZ_desir := RangeToMatrix(12,9,3,1);
      Q_desir := IK3(XYZ_desir);
      Qd1_desir := Mzeros(3,1);
      Qd2_desir := Mzeros(3,1);
  end;
  if InvDynON then begin

{******************************************************************************
 * INVERSE DYNAMICS - LABORATORIAL WORK
 * Implement here the inverse dynamics
 * Note that the variables required for the implementation are already created
 ******************************************************************************}

    // Computation of the matrices
    M_mat := Mzeros(3,3);
    C_mat := Mzeros(3,3);
    B_mat := Mzeros(3,3);
    if GetRCValue(10, 2) = 0 then begin
      Phi_mat := Mzeros(3,3);
    end else begin
      Phi_mat := Mzeros(3,3);
    end;

    // Debug
    MatrixToRangeF(24, 7, M_mat, '%.3f');
    MatrixToRangeF(28, 7, C_mat, '%.3f');
    MatrixToRangeF(32, 7, B_mat, '%.3f');
    MatrixToRangeF(24,10, Phi_mat, '%.3f');



    // Outter loop
    // - reference
    R_mat := Mzeros(3,1);
    // - output
    Aq_mat := Mzeros(3,1);

    // Debug
    MatrixToRangeF(16, 8, R_mat, '%.3f');
    MatrixToRangeF(16, 9, Aq_mat, '%.3f');



    // Inner loop - torque computation
    U_mat := Mzeros(3,1);
    
    // Debug
    MatrixToRangeF(16,10, U_mat, '%.3f');



    // Set torque of the joints
    SetTorque(U_mat);

    

    // Error computation
    Err_mat := MSub(Q_meas,Q_desir);

    // Debug
    SetRCValue(20, 8, format('%.3g',[Deg(MGetV(Err_mat,0,0))]));
    SetRCValue(21, 8, format('%.3g',[Deg(MGetV(Err_mat,1,0))]));
    SetRCValue(22, 8, format('%.3g',[MGetV(Err_mat,2,0)]));

{******************************************************************************
 ******************************************************************************}

  end;



  // SimTwo Sheet
  // - joint value (forward kinematics)
  SetRCValue(2, 8, format('%.3g',[Deg(GetAxisPos(0, iJ1Axis))]));
  SetRCValue(3, 8, format('%.3g',[Deg(GetAxisPos(0, iJ2Axis))]));
  SetRCValue(4, 8, format('%.3g',[GetAxisPos(0, iZAxis)]));

  // - pen position (horizontal/vertical SCARA)
  SetRCValue(2, 13, format('%.3g',[PosPen.x]));
  SetRCValue(12, 8, format('%.3g',[PosPen.x]));
  if GetRCValue(10, 2) = 0 then begin
    SetRCValue(3, 13, format('%.3g',[PosPen.y]));
    SetRCValue(4, 13, format('%.3g',[PosPen.z - 0.25]));
    SetRCValue(13, 8, format('%.3g',[PosPen.y]));
    SetRCValue(14, 8, format('%.3g',[PosPen.z - 0.25]));
  end else begin
    SetRCValue(3, 13, format('%.3g',[PosPen.z]));
    SetRCValue(4, 13, format('%.3g',[-PosPen.y - 0.25]));
    SetRCValue(13, 8, format('%.3g',[PosPen.z]));
    SetRCValue(14, 8, format('%.3g',[-PosPen.y - 0.25]));
  end;

  // - inverse dynamics activated
  if InvDynON then begin
    SetRCValue(10, 8, 'ON');
  end else begin
    SetRCValue(10, 8, 'OFF');
  end;
end;



{******************************************************************************
 * INITIALIZE
 * Initialization procedure is called only once when the script is started
 ******************************************************************************}
procedure Initialize;
begin
  iJ1Axis := GetAxisIndex(0, 'Joint1', 0);
  iJ2Axis := GetAxisIndex(0, 'Joint2', 0);
  iZAxis  := GetAxisIndex(0, 'SlideZ', 0);
  iPen := GetSolidIndex(0, 'Pen');

  Q_meas := Mzeros(3,1);
  Qd1_meas := Mzeros(3,1);
  Q_desir := Mzeros(3,1);
  Qd1_desir := Mzeros(3,1);
  Qd2_desir := Mzeros(3,1);
  W_mat := Mzeros(3,1);
  K0_mat := Mzeros(3,3);
  K1_mat := Mzeros(3,3);
end;
