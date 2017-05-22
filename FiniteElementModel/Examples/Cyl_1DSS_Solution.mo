within FiniteElementModel.Examples;
model Cyl_1DSS_Solution
  import Modelica.Math.*;
  import Modelica.SIunits.*;
  import Modelica.Constants.*;

parameter Length r_f = 0.005;
parameter Length r_g = 0.0055;
parameter Length r_c = 0.007;
parameter Length H_f = 1;

ThermalConductivity k_f;
ThermalConductivity k_g;
ThermalConductivity k_c;

ThermalResistance R_condg;
ThermalResistance R_condc;
ThermalResistance R_h;
 Temperature T_rf;
 Temperature T_SS(start = 500);

Real gppp(unit="W/m3") = Pin/(r_f^2*pi*H_f);

Length rr = 0;

  Modelica.Blocks.Interfaces.RealInput Pin(unit="W") "Inlet Power"
    annotation (Placement(transformation(extent={{-140,30},{-100,70}}),
        iconTransformation(extent={{-120,40},{-100,60}})));
  Modelica.Blocks.Interfaces.RealInput h(unit="W/(m2.K)")
    "Convection Coefficient"
    annotation (Placement(transformation(extent={{-140,-80},{-100,-40}}),
        iconTransformation(extent={{-120,-60},{-100,-40}})));
  Modelica.Blocks.Interfaces.RealInput T_inf(unit="K") "Coolant Temperature"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}}),
        iconTransformation(extent={{-120,-10},{-100,10}})));
  Modelica.Blocks.Interfaces.RealOutput Tmax(unit="K") "Centerline Temperature"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
algorithm

  (k_f,,) :=FiniteElementModel.FEM_Functions.Material_Fuel(T_SS);
  (k_g,,) :=FiniteElementModel.FEM_Functions.Material_Gap(T_SS);
  (k_c,,) :=FiniteElementModel.FEM_Functions.Material_Clad(T_SS);

  R_condg :=log(r_g/r_f)/(2*pi*k_g*H_f);
  R_condc :=log(r_c/r_g)/(2*pi*k_c*H_f);
  R_h :=1/(h*2*pi*r_c*H_f);

  T_rf :=gppp*pi*r_f^2*H_f*(R_condg + R_condc + R_h) + T_inf;

  T_SS :=gppp*r_f^2/(4*k_f)*(1 - (rr/r_f)^ 2) + T_rf;

  Tmax := max(T_SS);
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}})));
end Cyl_1DSS_Solution;
