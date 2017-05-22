within FiniteElementModel.FEM_Functions;
function Material_Gap

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.ThermalConductivity k;
  output Modelica.SIunits.Density rho;
  output Modelica.SIunits.SpecificHeatCapacity cp;

algorithm

  k :=0.1;
  rho :=200;
  cp :=4000;

end Material_Gap;
