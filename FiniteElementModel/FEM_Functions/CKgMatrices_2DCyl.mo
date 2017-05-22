within FiniteElementModel.FEM_Functions;
function CKgMatrices_2DCyl
  "Construct area integral matrices/vectors (C, K, and g) and volume scaled temperature"

  import Modelica.SIunits.*;
  import Modelica.Constants.*;
  import Modelica.Math.*;

  input Integer n_nodes;
  input Integer n_elem;
  input Integer n_subd;

  input Integer t[4,n_elem];
  input Real sd[1,n_subd];
  input Length p[2,n_nodes];
  input Temperature Temp[n_nodes];
  input Volume volFuel;
  input FiniteElementModel.FEM_Functions.Material_Fuel fuelType;
  input FiniteElementModel.FEM_Functions.Material_Gap gapType;
  input FiniteElementModel.FEM_Functions.Material_Clad claddingType;

  output HeatCapacity C[n_nodes,n_nodes];
  output ThermalConductance K[n_nodes,n_nodes];
  output Power g[n_nodes];

  output Temperature Tavg_f;

protected
  Integer i;
  Integer j;
  Integer k;
  Integer d;

  ThermalConductivity ke;
  Real g_pppe;
  Density rhoe;
  SpecificHeatCapacity cpe;

  Length ri, rj, rk;
  Length zi, zj, zk;

  Length rij, rik, rjk;
  Length zij, zik, zjk;
  Area bijk;

  Area Ae;

  Length cii, cij, cik;
  Length cji, cjj, cjk;
  Length cki, ckj, ckk;

  Length rce;
  Volume Ve;

  Area kii, kij, kik;
  Area kji, kjj, kjk;
  Area kki, kkj, kkk;

  Length gi, gj, gk;

  Length rr;
  Length zz;

  Real wi, wj, wk;
  Temperature Tempe;

algorithm
  Tavg_f :=0;

  for ie in 1:n_elem loop

    //%Specify nodes and subdomain associated with each element
    i :=t[1, ie];
    j :=t[2, ie];
    k :=t[3, ie];
    d :=t[4, ie];

    //%Specify properties of each element
    g_pppe :=sd[1, d];
//     ke     :=sd[1, d];
//     g_pppe :=sd[2, d];
//     rhoe   :=sd[3, d];
//     cpe    :=sd[4, d];

    //%Determine geometric characteristics of each element
    ri := p[1,i]; rj := p[1,j]; rk := p[1,k];
    zi := p[2,i]; zj := p[2,j]; zk := p[2,k];

    rij :=rj - ri; rik :=rk - ri; rjk :=rk - rj;
    zij :=zj - zi; zik :=zk - zi; zjk :=zk - zj;
    bijk :=rij*zjk - rjk*zij;
    assert(bijk>0,"bijk must be greater than 0");

    Ae :=abs(bijk)/2;

    cii := 6*ri + 2*rj + 2*rk; cij := 2*ri + 2*rj + rk;   cik := 2*ri + rj + 2*rk;
    cji := cij;                cjj := 2*ri + 6*rj + 2*rk; cjk := ri + 2*rj + 2*rk;
    cki := cik;                ckj := cjk;                ckk := 2*ri + 2*rj + 6*rk;

    rce := (ri + rj + rk)/3;

    Ve := 2*pi*rce*Ae;

    kii := rjk*rjk + zjk*zjk; kij := -(rjk*rik + zjk*zik); kik := rjk*rij + zjk*zij;
    kji := kij;               kjj := rik*rik + zik*zik;    kjk := -(rik*rij + zik*zij);
    kki := kik;               kkj := kjk;                  kkk := rij*rij + zij*zij;

    gi := 2*ri + rj + rk;     gj := ri + 2*rj + rk;      gk := ri + rj + 2*rk;

    rr :=(ri + rj + rk)/3;
    zz :=(zi + zj + zk)/3;

    wi :=1/bijk*((rj*zk - rk*zj) - rr*zjk + rjk*zz);
    wj :=1/bijk*((rk*zi - ri*zk) + rr*zik - rik*zz);
    wk :=1/bijk*((ri*zj - rj*zi) - rr*zij + rij*zz);

    Tempe :=wi*Temp[i] + wj*Temp[j] + wk*Temp[k];

    //Temperature dependent materials properties
    if d == n_subd then
      //Cladding
      (ke,rhoe,cpe) :=claddingType(Tempe);

    elseif d == n_subd-1 then
      //Gap
      (ke,rhoe,cpe) :=gapType(Tempe);
    else
      //Fuel
      (ke,rhoe,cpe) :=fuelType(Tempe);
      Tavg_f :=Tavg_f + Tempe*Ve/volFuel;
    end if;

    //%Construct capacitance matrix (C)
    C[i,i] :=C[i, i] + cii*2*pi*rhoe*cpe*Ae/60;
    C[i,j] :=C[i, j] + cij*2*pi*rhoe*cpe*Ae/60;
    C[i,k] :=C[i, k] + cik*2*pi*rhoe*cpe*Ae/60;
    C[j,i] :=C[j, i] + cji*2*pi*rhoe*cpe*Ae/60;
    C[j,j] :=C[j, j] + cjj*2*pi*rhoe*cpe*Ae/60;
    C[j,k] :=C[j, k] + cjk*2*pi*rhoe*cpe*Ae/60;
    C[k,i] :=C[k, i] + cki*2*pi*rhoe*cpe*Ae/60;
    C[k,j] :=C[k, j] + ckj*2*pi*rhoe*cpe*Ae/60;
    C[k,k] :=C[k, k] + ckk*2*pi*rhoe*cpe*Ae/60;

    //%Construct conductance matrix (K)
    K[i,i] :=K[i, i] + kii*2*pi*rce*ke/(4*Ae);
    K[i,j] :=K[i, j] + kij*2*pi*rce*ke/(4*Ae);
    K[i,k] :=K[i, k] + kik*2*pi*rce*ke/(4*Ae);
    K[j,i] :=K[j, i] + kji*2*pi*rce*ke/(4*Ae);
    K[j,j] :=K[j, j] + kjj*2*pi*rce*ke/(4*Ae);
    K[j,k] :=K[j, k] + kjk*2*pi*rce*ke/(4*Ae);
    K[k,i] :=K[k, i] + kki*2*pi*rce*ke/(4*Ae);
    K[k,j] :=K[k, j] + kkj*2*pi*rce*ke/(4*Ae);
    K[k,k] :=K[k, k] + kkk*2*pi*rce*ke/(4*Ae);

    //%Construct generation vector (g)
    g[i] :=g[i] + gi*2*pi*g_pppe*Ae/12;
    g[j] :=g[j] + gj*2*pi*g_pppe*Ae/12;
    g[k] :=g[k] + gk*2*pi*g_pppe*Ae/12;
  end for;

end CKgMatrices_2DCyl;
