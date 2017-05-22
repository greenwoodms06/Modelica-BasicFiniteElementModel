within FiniteElementModel.FEM_Functions;
function Material_Fuel

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.ThermalConductivity k;
  output Modelica.SIunits.Density rho;
  output Modelica.SIunits.SpecificHeatCapacity cp;

algorithm

  k :=1;
  rho :=1000;
  cp :=4000;

end Material_Fuel;
