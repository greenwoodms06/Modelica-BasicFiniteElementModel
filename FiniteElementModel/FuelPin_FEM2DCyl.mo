within FiniteElementModel;
model FuelPin_FEM2DCyl
  import FiniteElementModel;
  //Finite Element solution of a 2D Cylinder with three distinct radial zones

import Modelica.SIunits.*;
import Modelica.Constants.*;
import Modelica.Math.*;
import Modelica.Fluid.Types.Dynamics;

  // Material Type Selection
replaceable function fuelType =
      FiniteElementModel.FEM_Functions.Material_Fuel
                                  constrainedby
    FiniteElementModel.FEM_Functions.Material_Fuel
   annotation(choicesAllMatching=true,Dialog(group="Materials"));

replaceable function gapType =
      FiniteElementModel.FEM_Functions.Material_Gap
                                 constrainedby
    FiniteElementModel.FEM_Functions.Material_Gap
   annotation(choicesAllMatching=true,Dialog(group="Materials"));

replaceable function claddingType =
      FiniteElementModel.FEM_Functions.Material_Clad
                                  constrainedby
    FiniteElementModel.FEM_Functions.Material_Clad
   annotation(choicesAllMatching=true,Dialog(group="Materials"));

  // Fuel Pin Geometry
  parameter Length r_0f = 0 "Inner radius of fuel"
    annotation (Dialog(tab="General", group="Geometry"));
  parameter Length r_f = 0.005 "Outer radius of fuel element"
    annotation (Dialog(tab="General", group="Geometry"));
  parameter Length r_g = 0.0055 "Outer radius of gap"
    annotation (Dialog(tab="General", group="Geometry"));
  parameter Length r_c = 0.007 "Outer radius of cladding"
    annotation (Dialog(tab="General", group="Geometry"));
  parameter Length H_f = 1 "Height of fuel element"
  annotation (Dialog(tab="General", group="Geometry"));

  // Nodalization Scheme
  parameter Integer n_finner = 2 "Inner # of primary radial nodes in fuel"
    annotation (Dialog(tab="General", group="Nodalization"));
  parameter Integer n_ginner = 1 "Inner # of primary radial nodes in gap"
    annotation (Dialog(tab="General", group="Nodalization"));
  parameter Integer n_cinner = 1 "Inner # of primary radial nodes in cladding"
    annotation (Dialog(tab="General", group="Nodalization"));
  parameter Integer nNodes = 3 "# of discrete axial volume nodes"
    annotation (Dialog(tab="General", group="Nodalization"));

  // Total # of fuel pins for determination of total power output
  parameter Integer nFuelPins = 1 "Total # of fuel pin elements"
    annotation (Dialog(tab="General",group="Power Scaling"));

  // Determine # of nodes
  final parameter Integer n_zinner = nNodes - 1
    "Inner # of primary nodes in axial direction";
  final parameter Integer n_radial = 1 + n_finner + 1 + n_ginner + 1 + n_cinner + 1
    "Total # of nodes in radial direction";
  final parameter Integer n_axial = 1 + n_zinner + 1
    "Total # of primary nodes in axial direction";
  final parameter Integer n_nodes = n_radial*n_axial + (n_radial-1)*(n_axial-1)
    "# of nodes";
  final parameter Integer n_elem = 4*(n_radial-1)*(n_axial-1) "# of elements";
  final parameter Integer n_edges = 2*(n_radial-1) + 2*(n_axial-1) "# of edges";

  final parameter Integer n_subd = (n_axial-1) + 2 "# of subdomains";
  final parameter Integer n_bound = 1 + (n_axial-1) + 2 "# of boundaries";

  // Determine required volumes
  final parameter Volume volFuel = H_f*pi*(r_f^2-r_0f^2);
  final parameter Volume volFuelNode = volFuel/(n_axial-1)
    "Volume of fuel in an axial nodal plane";

  // Initial condition parameters
  parameter Temperature T_0fuel = 772 "Initial Fuel Temperature"
    annotation (Dialog(tab="Initial Conditions", group="Temperature"));
  parameter Temperature T_0gap = 722 "Initial Gap Temperature"
    annotation (Dialog(tab="Initial Conditions", group="Temperature"));
  parameter Temperature T_0clad = 510 "Initial Cladding Temperature"
    annotation (Dialog(tab="Initial Conditions", group="Temperature"));

  Temperature[n_axial*(n_finner+2)] T_0f = T_0fuel*ones(n_axial*(n_finner+2));
  Temperature[n_axial*(n_ginner+1)] T_0g = T_0gap*ones(n_axial*(n_ginner+1));
  Temperature[n_axial*(n_cinner+1)] T_0c = T_0clad*ones(n_axial*(n_cinner+1));

  // Solution matrix (Temp) initialization
  Temperature Temp[n_nodes] "Solution";

  // Intermediate parameters of interest
  Temperature Tmax "Maximum Temperature in fuel element";
  Temperature T_wall[n_axial] "Cladding outer wall temperatures";
  Power Q_flow[n_axial-1] "Heat flow out at cladding/coolant (c/c) interface";
  Power Q_total "Total heat flow due to convection at c/c interface";
  Temperature T_eff_f "Fuel effective (average) temperature in the fuel region";

  /*=============================*/
  /* Initial Condition Option    */
  /*=============================*/

  parameter Dynamics energyDynamics = Dynamics.DynamicFreeInitial
    "Formulation of energy balances";

  /* =================== */
  /* Interactive Options */
  /* =================== */

  // Inputs
  parameter Boolean use_Power_in = false
    "Check to use the axial power distribution from the point kinetics model"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true));

  Modelica.Blocks.Interfaces.RealInput[n_axial-1] Power_in(unit="W") if use_Power_in
    "Axial power distribution from the point kinetics model"
    annotation (Placement(transformation(extent={{-19.5,-19.5},{19.5,19.5}},
          rotation=0,
        origin={-119.5,30.5}),
                       iconTransformation(extent={{7.375,-7.875},{-7.375,7.875}},
        rotation=180,
        origin={-17.375,29.875})));

  parameter Boolean use_h_conv_in = false
    "Check to use heat convection coefficient from the pipe model"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true));

  Modelica.Blocks.Interfaces.RealInput h_conv_in[n_axial-1](unit="W/(m2.K)") if use_h_conv_in
    "Input heat convection coefficient from the pipe model"
    annotation (Placement(transformation(extent={{-140,-50},{-100,-10}},
          rotation=0), iconTransformation(extent={{7.25,-7.25},{-7.25,7.25}},
        rotation=180,
        origin={-17.75,-17.75})));

  parameter Boolean use_T_coolant_in = false
    "Check to use the coolant temperature from the pipe model"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true));

  Modelica.Blocks.Interfaces.RealInput T_coolant_in[n_axial-1](unit="K") if use_T_coolant_in
    "Input the coolant temperature from the pipe model"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
          rotation=0), iconTransformation(extent={{7.25,-7.25},{-7.25,7.25}},
        rotation=180,
        origin={-18.25,4.75})));

  // Outputs
  parameter Boolean use_Teff_out = false
    "Check to send fuel effective temperature to point kinetics model"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true));

  Modelica.Blocks.Interfaces.RealOutput Teff_out(unit="K") if use_Teff_out
    "Output the fuel effective temperature to point kinetics model"
    annotation (Placement(transformation(extent={{100,-40},{140,0}},
          rotation=0), iconTransformation(extent={{7.25,-7.25},{-7.25,7.25}},
        rotation=180,
        origin={17.75,-10.25})));

  parameter Boolean use_Power_out = false
    "Check to send the heat flow to the pipe model"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true));

  Modelica.Blocks.Interfaces.RealOutput[n_axial-1] Power_out(unit="W") if use_Power_out
    "Output the heat flow to the pipe model"
    annotation (Placement(transformation(extent={{-19.5,-19.5},{19.5,19.5}},
          rotation=0,
        origin={119.5,20.5}),
                       iconTransformation(extent={{7.5,-7.75},{-7.5,7.75}},
        rotation=180,
        origin={18.5,19.75})));

protected
  Modelica.Blocks.Interfaces.RealInput[n_axial-1] Pin_Internal
    "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.RealInput[n_axial-1] h_conv_Internal(unit="W/(m2.K)")
    "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.RealInput[n_axial-1] T_coolant_Internal(unit="K")
    "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.RealOutput Teff_Internal(unit="K")
    "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.RealOutput[n_axial-1] Pout_Internal
    "Needed to connect to conditional connector";

  // Point (p), edge (e), and triangle (t) matrix initialization
  Length p[2,n_nodes];
  Integer e[3,n_edges];
  Integer t[4,n_elem];

  // Subdomain (sd) matrix initialization
  //  k, rho, and cp moved to CKgMatrices function
  //  ThermalConductivity k[1,n_subd];
  //  Density rho[1,n_subd];
  //  HeatCapacity cp[1,n_subd];
  Real g_pppf[n_axial-1](unit="W/m3");
  Real g_ppp[1,n_subd](unit="W/m3");
  Real sd[1,n_subd];

  // Boundary condition (bc) matrix initialization
  CoefficientOfHeatTransfer h[1,n_bound];
  HeatFlux q_pp_s[1,n_bound];
  Temperature T_infinity[1,n_bound];
  Real bc[3,n_bound];

  // Boundary conditions matrix initialization (H, qs, and Htinf)
  ThermalConductance H[n_nodes,n_nodes];
  Power qs[n_nodes];
  Power Htinf[n_nodes];

  // Capacitance (C), conduction (K), and generation (g) matrix initialization
  HeatCapacity C[n_nodes,n_nodes];
  ThermalConductance K[n_nodes,n_nodes];
  Power g[n_nodes];

  Temperature[n_nodes] Temp_0;

algorithm
  (p,e,t,Temp_0) :=FiniteElementModel.FEM_Functions.petMatrices_2DCyl(
      r_0f=r_0f,
      r_f=r_f,
      r_g=r_g,
      r_c=r_c,
      H_f=H_f,
      n_finner=n_finner,
      n_ginner=n_ginner,
      n_cinner=n_cinner,
      n_zinner=n_zinner,
      n_radial=n_radial,
      n_axial=n_axial,
      n_nodes=n_nodes,
      n_elem=n_elem,
      n_edges=n_edges,
      n_subd=n_subd,
      n_bound=n_bound,
      T_0f=T_0f,
      T_0g=T_0g,
      T_0c=T_0c);

initial equation
  if energyDynamics == Dynamics.SteadyStateInitial then
    der(Temp)   = zeros(n_nodes);
  else
    Temp = Temp_0;
  end if;

equation

  /* =================== */
  /* Interactive Options */
  /* =================== */
  /* 
  Connect input/output variables to the internally used versions of the same
  variable. If external input variables are not used, default internal values
  as specified are used.
  */

  connect(Pout_Internal, Power_out);
  Pout_Internal = Q_flow;

  // Extract wall temperature values
  T_wall = {Temp[n_radial+(2*n_radial-1)*(i-1)] for i in 1:n_axial};

  // Calculate heat flux at each wall edge based on convection and fluid temp
  Q_flow[1:n_axial-1] = nFuelPins*{pi*2*r_c*H_f/(n_axial-1)*h_conv_Internal[i]*((T_wall[i]+T_wall[i+1])/2 - T_coolant_Internal[i]) for i in 1:n_axial-1};
  Q_total = sum(Q_flow);

  connect(Teff_Internal, Teff_out);
  Teff_Internal = T_eff_f;

  // ***************************************************************************

  /*
  Construct subdomain matrix (sd) 4 x n_subd. Each matrix (i.e. k, g_ppp, 
  rho, and cp) properties are defined as [fuel,gap,cladding].
  */

  // Axial energy generation (W) in fuel
  connect(Pin_Internal, Power_in);
  if not use_Power_in then
    // Uniform power profile with a total power/pin as stated
//      Pin_Internal = zeros(n_axial-1);
      Pin_Internal = (10000/(n_axial-1))*ones(n_axial-1);
  end if;
  g_pppf = Pin_Internal/volFuelNode;
  g_ppp = [{g_pppf},0,0];
  sd = [g_ppp];

  // properties k, rho, and cp for Fuel, Gap, and Cladding are located in
  // CKgMatrices_2DCyl as they are temperature dependent. Therefore sd only
  // contains g_ppp.

// k =     [{2.3*ones(n_axial-1)},0.1513,15];
// g_ppp = [g_pppf,0,0];
// rho =   [{11000*ones(n_axial-1)},0.164,8000];
// cp =    [{300*ones(n_axial-1)},5193,466];
//
//   sd = [k;g_ppp;rho;cp];

  // ***************************************************************************

  /*
  Construct boundary matrix (bc) 3 x n_bound. Each matrix (i.e. h, q_pp_s, 
  T_infinity) associates with a boundary in the following order 
  [bottom,cladding,top,centerline].
  */

  //Convective boundaries: Convection Heat Transfer Coefficient
  connect(h_conv_Internal, h_conv_in);
  if not use_h_conv_in then
//      h_conv_Internal = zeros(n_axial-1);
      h_conv_Internal = 2e4 * ones(n_axial-1);
  end if;
  h = [0,{h_conv_Internal},0,0];

  //Convective boundaries: Coolant Temperature
  connect(T_coolant_Internal, T_coolant_in);
  if not use_T_coolant_in then
    T_coolant_Internal = 500*ones(n_axial-1);
  end if;
  T_infinity = [0,{T_coolant_Internal},0,0];

  //Specified Heat Flux boundary (can be + or -)
  q_pp_s = [0,zeros(1,n_axial-1),0,0];

  //Boundary matrix
  bc = [h;T_infinity;q_pp_s];

  // ***************************************************************************

  /* =============== */
  /* FEMethod Solver */
  /* =============== */

  (C,K,g,T_eff_f) =FiniteElementModel.FEM_Functions.CKgMatrices_2DCyl(
      n_nodes=n_nodes,
      n_elem=n_elem,
      n_subd=n_subd,
      t=t,
      sd=sd,
      p=p,
      Temp=Temp,
      volFuel=volFuel,
      fuelType=fuelType,
      gapType=gapType,
      claddingType=claddingType);

  (H,qs,Htinf) =FiniteElementModel.FEM_Functions.BoundaryMatrices_2DCyl(
      n_nodes=n_nodes,
      n_edges=n_edges,
      n_bound=n_bound,
      e=e,
      p=p,
      bc=bc);

  // Time Dependent Solution
   C*der(Temp) + (K + H)*Temp = g + qs + Htinf;

  // Steady State Solution
  // S*Temp = g + qs + Htinf;

  Tmax = max(Temp);
  annotation (Documentation(info="<html>
<p><br>Finite element (Galerkin weighted residual method) solution to a 2D cylinder with user specified properties, internal heat generation, and external heat transfer.</p>
<p><br>The following figure represents a basic guideline defining what occurs in the solver.</p>
</html>"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics={
        Rectangle(
          extent={{-50,60},{50,-60}},
          lineColor={113,113,113},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={185,185,185}),
        Rectangle(
          extent={{-40,60},{40,-60}},
          lineColor={170,255,85},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={255,170,85}),
        Rectangle(
          extent={{-32,60},{34,-60}},
          lineColor={63,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={216,0,0}),
        Text(
          extent={{-12,-68},{12,-78}},
          lineColor={63,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={216,0,0},
          textString="Fuel",
          fontSize=16),
        Text(
          extent={{20,-68},{44,-78}},
          lineColor={63,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={216,0,0},
          fontSize=16,
          textString="Gap"),
        Text(
          extent={{42,-68},{66,-78}},
          lineColor={63,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={216,0,0},
          fontSize=16,
          textString="Cladding"),
        Line(points={{0,-70},{0,-64},{2,-64},{0,-60},{-2,-64},{0,-64}}, color={63,
              0,0}),
        Line(points={{36,-70},{36,-64},{38,-64},{36,-60},{34,-64},{36,-64}},
            color={63,0,0}),
        Line(points={{48,-70},{48,-64},{50,-64},{48,-60},{46,-64},{48,-64}},
            color={63,0,0}),
        Rectangle(
          extent={{50,60},{64,-60}},
          fillPattern=FillPattern.Solid,
          fillColor={0,128,255},
          pattern=LinePattern.None),
        Rectangle(
          extent={{-64,60},{-50,-60}},
          fillPattern=FillPattern.Solid,
          fillColor={0,128,255},
          pattern=LinePattern.None),
        Text(
          extent={{74,-50},{98,-60}},
          lineColor={63,0,0},
          fillPattern=FillPattern.VerticalCylinder,
          fillColor={216,0,0},
          fontSize=16,
          textString="Coolant"),
        Line(
          points={{0,-5},{0,1},{2,1},{0,5},{-2,1},{0,1}},
          color={63,0,0},
          origin={72,-55},
          rotation=90)}),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),
        graphics={
        Text(
          visible=use_Teff_out,
          extent={{0,-7.5},{26,7.5}},
          lineColor={255,0,0},
          origin={22,-1.5},
          rotation=0,
          textString="Teff"),
        Text(
          visible=use_h_conv_in,
          extent={{0,-6},{28,6}},
          lineColor={0,0,255},
          origin={-54,-10},
          rotation=0,
          textString="hconv"),
        Text(
          visible=use_Power_in,
          extent={{0,-6.5},{20,6.5}},
          lineColor={0,0,255},
          origin={-46,41.5},
          rotation=0,
          textString="Pin"),
        Text(
          visible=use_T_coolant_in,
          extent={{0,-10},{31,10}},
          lineColor={0,0,255},
          origin={-60,14},
          rotation=0,
          textString="Tcoolant"),
        Text(
          visible=use_Power_out,
          extent={{0,-7},{24,7}},
          lineColor={255,0,0},
          origin={26,31},
          rotation=0,
          textString="Pout"),
        Text(
          extent={{6,-82},{-4,-54}},
          lineColor={0,0,0},
          textString="%nNodes"),
        Text(
          extent={{-24,106},{24,92}},
          lineColor={28,108,200},
          textString="%name")}));
end FuelPin_FEM2DCyl;
