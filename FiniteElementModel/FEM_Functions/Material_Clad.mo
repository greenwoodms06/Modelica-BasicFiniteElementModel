within FiniteElementModel.FEM_Functions;
function Material_Clad

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.ThermalConductivity k;
  output Modelica.SIunits.Density rho;
  output Modelica.SIunits.SpecificHeatCapacity cp;

algorithm

  k :=10;
  rho :=1000;
  cp :=4000;

end Material_Clad;
