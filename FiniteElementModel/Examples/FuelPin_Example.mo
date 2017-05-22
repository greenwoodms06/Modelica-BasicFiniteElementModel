within FiniteElementModel.Examples;
model FuelPin_Example

  Cyl_1DSS_Solution cyl_1DSS_Solution
    annotation (Placement(transformation(extent={{-16,-40},{4,-20}})));
  FuelPin_FEM2DCyl fuelPin_FEM2DCyl(
    use_Power_in=true,
    use_h_conv_in=true,
    use_T_coolant_in=true,
    use_Teff_out=true,
    use_Power_out=true,
    redeclare function fuelType = FEM_Functions.Material_Fuel,
    redeclare function gapType = FEM_Functions.Material_Gap,
    redeclare function claddingType = FEM_Functions.Material_Clad)
    annotation (Placement(transformation(extent={{-36,0},{28,58}})));

  Modelica.Blocks.Sources.Constant P[3](k=10000/3*ones(3))
    annotation (Placement(transformation(extent={{-76,44},{-64,56}})));
  Modelica.Blocks.Sources.Constant T[3](k=500*ones(3))
    annotation (Placement(transformation(extent={{-76,24},{-64,36}})));
  Modelica.Blocks.Sources.Constant h[3](k=2e4*ones(3))
    annotation (Placement(transformation(extent={{-76,4},{-64,16}})));
  Modelica.Blocks.Sources.Constant P1(k=10000)
    annotation (Placement(transformation(extent={{-78,-31},{-66,-19}})));
equation
  connect(h.y, fuelPin_FEM2DCyl.h_conv_in) annotation (Line(points={{-63.4,10},
          {-34,10},{-34,23.8525},{-9.68,23.8525}}, color={0,0,127}));
  connect(T.y, fuelPin_FEM2DCyl.T_coolant_in) annotation (Line(points={{-63.4,
          30},{-9.84,30},{-9.84,30.3775}}, color={0,0,127}));
  connect(P.y, fuelPin_FEM2DCyl.Power_in) annotation (Line(points={{-63.4,
          50},{-34,50},{-34,37.6637},{-9.56,37.6637}},
                                                  color={0,0,127}));
  connect(cyl_1DSS_Solution.T_inf, T[1].y) annotation (Line(points={{-17,-30},{-50,
          -30},{-50,30},{-63.4,30}}, color={0,0,127}));
  connect(cyl_1DSS_Solution.h, h[1].y) annotation (Line(points={{-17,-35},{-54,-35},
          {-54,10},{-63.4,10}}, color={0,0,127}));
  connect(P1.y, cyl_1DSS_Solution.Pin) annotation (Line(points={{-65.4,-25},{-42,
          -25},{-17,-25}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})), Documentation(info="<html>
<p>This example shows that the 1D centerline temperature matches the Finite Element solution with the same inputs and constant thermal conductivity properties.</p>
<p>To use make sure the dimensions of the two models are the same and go the the materials property functions and make them constant and identical to the 1D model.</p>
</html>"),
    experiment(StopTime=500),
    __Dymola_experimentSetupOutput);
end FuelPin_Example;
