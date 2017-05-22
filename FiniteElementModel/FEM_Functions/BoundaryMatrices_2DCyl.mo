within FiniteElementModel.FEM_Functions;
function BoundaryMatrices_2DCyl
  "Construct bondary integral matrices (qs, H, Htinf)"

  import Modelica.SIunits.*;
  import Modelica.Constants.*;

  input Integer n_nodes;
  input Integer n_edges;
  input Integer n_bound;

  input Integer e[3,n_edges];
  input Length p[2,n_nodes];
  input Real bc[3,n_bound];

  output ThermalConductance H[n_nodes,n_nodes];
  output Power qs[n_nodes];
  output Power Htinf[n_nodes];

protected
  Integer i;
  Integer j;
  Integer bnd;

  Length ri, rj;
  Length zi, zj;

  CoefficientOfHeatTransfer hb;
  Temperature Tinfb;
  HeatFlux qs_ppb;

  Length rij;
  Length zij;
  Length sij;

  Length qi, qj;
  Length hi, hj;
  Length hii, hij, hji, hjj;

algorithm
for ib in 1:n_edges loop

    i := e[1, ib];
    j := e[2, ib];
    bnd := e[3, ib];

    ri := p[1, i];  rj := p[1, j];
    zi := p[2, i];  zj := p[2, j];

    hb     := bc[1, bnd];
    Tinfb  := bc[2, bnd];
    qs_ppb := bc[3, bnd];

    rij := rj - ri;
    zij := zj - zi;
    sij := sqrt(rij^2 + zij^2);

    qi := 2*ri + rj;  qj := ri + 2*rj;

    hi := 2*ri + rj;  hj := ri + 2*rj;

    hii := 3*ri + rj;  hij := ri + rj;  hji := hij; hjj := ri + 3*rj;

    //%Construct specified flux vector (qs)
    qs[i] :=qs[i] + qi*2*pi*qs_ppb*sij/6;
    qs[j] :=qs[j] + qj*2*pi*qs_ppb*sij/6;

    //%Construct convective matrix (H)
    H[i,i] :=H[i, i] + hii*2*pi*hb*sij/12;
    H[i,j] :=H[i, j] + hij*2*pi*hb*sij/12;
    H[j,i] :=H[j, i] + hji*2*pi*hb*sij/12;
    H[j,j] :=H[j, j] + hjj*2*pi*hb*sij/12;

    //%Construct Tinfinity convective matrix (Htinf)
    Htinf[i] :=Htinf[i] + hi*2*pi*hb*Tinfb*sij/6;
    Htinf[j] :=Htinf[j] + hj*2*pi*hb*Tinfb*sij/6;

end for;

end BoundaryMatrices_2DCyl;
