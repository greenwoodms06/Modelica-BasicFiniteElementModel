within FiniteElementModel.FEM_Functions;
function petMatrices_2DCyl
  "Mesh of a fuel element with a gas gap and cladding. A node is added between 'main' nodes."

  import Modelica.SIunits.*;

  input Length r_0f "Inner radius of fuel";
  input Length r_f "Outer radius of fuel element";
  input Length r_g "Outer radius of gap";
  input Length r_c "Outer radius of cladding";
  input Length H_f "Height of fuel element";

  input Integer n_finner "Inner # of primary radial nodes in fuel";
  input Integer n_ginner "Inner # of primary radial nodes in gap";
  input Integer n_cinner "Inner # of primary radial nodes in cladding";
  input Integer n_zinner "Inner # of primary nodes in axial direction";

  input Integer n_radial "Total # of nodes in radial direction";
  input Integer n_axial "Total # of primary nodes in axial direction";
  input Integer n_nodes "# of nodes";
  input Integer n_elem "# of elements";
  input Integer n_edges "# of edges";

  input Integer n_subd "# of subdomains";
  input Integer n_bound "# of boundaries";

  input Temperature[n_axial*(n_finner+2)] T_0f;
  input Temperature[n_axial*(n_ginner+1)] T_0g;
  input Temperature[n_axial*(n_cinner+1)] T_0c;

  output Length p[2,n_nodes];
  output Integer e[3,n_edges];
  output Integer t[4,n_elem];
  output Temperature[n_nodes] Temp_0;

protected
  Length dr_f = (r_f-r_0f)/(n_finner+1) "Radial spacing in fuel";
  Length dr_g = (r_g-r_f)/(n_ginner+1) "Radial spacing in gap";
  Length dr_c = (r_c-r_g)/(n_cinner+1) "Radial spacing in clad";
  Length dz = H_f/(n_zinner+1) "Axial spacing in fuel element";

  Length r_position;
  Length z_position;

  Length pr[1,n_nodes] "r location of node";
  Length pz[1,n_nodes] "z location of node";

  Integer ei[1,n_edges] "Left edge node";
  Integer ej[1,n_edges] "Right edge node";
  Integer ee[1,n_edges] "Boundary #";

  // Construct the triangle matrix (t): 4 x n_elem
  Integer ti[1,n_elem]
    "On Up/Down Triangles: most left node; On Left/Right Triangles: top node";
  Integer tj[1,n_elem]
    "First point counter-clockwise from point specified in ti";
  Integer tk[1,n_elem]
    "Second point counter-clockwise from point specified in ti";
  Integer td[1,n_elem] "Subdomain";

  Integer k_sd;
  Integer t_subd;
  Integer jshift;

  Temperature Tempn;

algorithm
/*--------------------------------------------------*/

  // Construct the point matrix (p) 2 x n_nodes
  // primary nodes
  r_position :=0;
  z_position :=0;
  for j in 1:n_axial loop
    for i in 1:n_radial loop

      if i <= n_finner + 2 then
        r_position :=r_0f + dr_f*(i - 1);
        Tempn :=T_0f[i+(j-1)*(n_finner+2)];

      elseif (i <= n_finner + n_ginner + 3 and i > n_finner + 2) then
        r_position :=r_f + dr_g*(i - (n_finner + 2));
        Tempn :=T_0g[(i-n_finner-2)+(j-1)*(n_ginner+1)];

      else
        r_position :=r_g + dr_c*(i - (n_finner + n_ginner + 3));
        Tempn :=T_0c[(i-n_finner-n_ginner-3)+(j-1)*(n_cinner+1)];

      end if;

      z_position :=dz*(j - 1);

      pr[1,i + (2*n_radial-1)*(j-1)] :=r_position;
      pz[1,i + (2*n_radial-1)*(j-1)] :=z_position;

      Temp_0[i+(2*n_radial-1)*(j-1)] :=Tempn;
    end for;
  end for;

  // Interstitial nodes
  r_position :=0;
  z_position :=0;
  for j in 1:n_axial-1 loop
    for i in 1:n_radial-1 loop

      if i <= n_finner + 1 then

        r_position :=r_0f + dr_f*(i - 1) + dr_f/2;
        Tempn :=T_0f[i+(j-1)*(n_finner+2)];

      elseif (i <= n_finner + n_ginner + 2 and i > n_finner + 1) then

        r_position :=r_f - dr_g/2 + dr_g*(i - (n_finner + 1));

        if i == n_finner + 2 then
          Tempn :=(T_0g[(i-n_finner-1)+(j-1)*(n_ginner+1)] + T_0f[i+(j-1)*(n_finner+2)])/2;
        else
          Tempn :=T_0g[(i-n_finner-1)+(j-1)*(n_ginner+1)];
        end if;

      else

        r_position :=r_g - dr_c/2 + dr_c*(i - (n_finner + n_ginner + 2));

        if i ==n_finner + n_ginner + 3 then
          Tempn :=(T_0c[(i-n_finner-n_ginner-2)+(j-1)*(n_cinner+1)] + T_0g[(i-n_finner-1)+(j-1)*(n_ginner+1)])/2;
        else
          Tempn :=T_0c[(i-n_finner-n_ginner-2)+(j-1)*(n_cinner+1)];
        end if;

      end if;

      z_position :=dz*(j - 1) + dz/2;

      pr[1,i + n_radial + (2*n_radial-1)*(j-1)] :=r_position;
      pz[1,i + n_radial + (2*n_radial-1)*(j-1)] :=z_position;

      Temp_0[i + n_radial + (2*n_radial-1)*(j-1)] :=Tempn;
    end for;
  end for;
  p :=[pr; pz];

/*--------------------------------------------------*/

   // Construct the edge matrix (e): 3 x n_edges
  for i in 1:n_edges loop
    if i<=n_radial-1 then
      // Bottom edge
      ei[1,i] :=i;
      ej[1,i] :=i + 1;
      ee[1,i] :=1;
    elseif i<=n_radial+n_axial-2 then
      // Right edge
      ei[1,i] :=i + (i - n_radial)*(2*n_radial - 1) - (i - n_radial);
      ej[1,i] :=i + (i + 1 - n_radial)*(2*n_radial - 1) - (i - n_radial);
      ee[1,i] :=ee[1,i-1]+1;
    elseif i<=2*n_radial+n_axial-3 then
      ei[1,i] :=n_nodes - (i - (n_radial + n_axial - 1));
      ej[1,i] :=n_nodes - (i - (n_radial + n_axial - 2));
      ee[1,i] :=n_axial + 1;
    else
      ei[1,i] :=n_nodes - (i - (n_radial + n_axial - 1)) + (i - (2*n_radial +
        n_axial - 2))*(1 - (2*n_radial - 1));
      ej[1,i] :=n_nodes - (i + 1 - (n_radial + n_axial - 1)) + (i + 1 - (2*
        n_radial + n_axial - 2))*(1 - (2*n_radial - 1));
      ee[1,i] :=n_axial + 1;
    end if;
  end for;
  e :=[ei; ej; ee];

/*--------------------------------------------------*/

  // Construct the triangle matrix (t): 4 x n_elem
  k_sd :=0;
  for j in 1:n_axial-1 loop
    for i in 1:n_radial-1 loop
      // Determine subdomain of element
      if i <= n_finner + 1 then
        t_subd :=k_sd + 1;
      elseif (i <= n_finner + n_ginner + 2 and i > n_finner + 1) then
        t_subd :=n_subd-1;
      else
        t_subd :=n_subd;
      end if;

      jshift :=4*(j - 1)*(n_radial - 1);
      // Up facing triangles
      ti[1,i + jshift] :=i + (j - 1)*(2*n_radial - 1);
      tj[1,i + jshift] :=i + 1 + (j - 1)*(2*n_radial - 1);
      tk[1,i + jshift] :=i + n_radial + (j - 1)*(2*n_radial - 1);
      td[1,i + jshift] :=t_subd;

      // Right facing triangles
      ti[1,n_radial + 2*(i-1) + jshift] :=i + (2*n_radial - 1) + (j - 1)*(2*
        n_radial - 1);
      tj[1,n_radial + 2*(i-1) + jshift] :=i + (j - 1)*(2*n_radial - 1);
      tk[1,n_radial + 2*(i-1) + jshift] :=i + n_radial + (j - 1)*(2*n_radial - 1);
      td[1,n_radial + 2*(i-1) + jshift] :=t_subd;

      // Left facing triangles
      ti[1,n_radial - 1 + 2*i + jshift] :=i + 2*n_radial + (j - 1)*(2*n_radial -
        1);
      tj[1,n_radial - 1 + 2*i + jshift] :=i + n_radial + (j - 1)*(2*n_radial - 1);
      tk[1,n_radial - 1 + 2*i + jshift] :=i + 1 + (j - 1)*(2*n_radial - 1);
      td[1,n_radial - 1 + 2*i + jshift] :=t_subd;

      // Down facing triangles
      ti[1,3*n_radial + (i-3) + jshift] :=i + (2*n_radial - 1) + (j - 1)*(2*
        n_radial - 1);
      tj[1,3*n_radial + (i-3) + jshift] :=i + n_radial + (j - 1)*(2*n_radial - 1);
      tk[1,3*n_radial + (i-3) + jshift] :=i + 2*n_radial + (j - 1)*(2*n_radial -
        1);
      td[1,3*n_radial + (i-3) + jshift] :=t_subd;

    end for;
  end for;
  t :=[ti; tj; tk; td];

/*--------------------------------------------------*/

end petMatrices_2DCyl;
