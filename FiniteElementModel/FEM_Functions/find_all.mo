within FiniteElementModel.FEM_Functions;
function find_all "Find all elements (e) in a vector"
  extends Modelica.Icons.Function;
  input Real e "Search for e";
  input Real v[:] "Integer vector";
  input Real eps(min=0) = 0
    "Element e is equal to a element v[i] of vector v if abs(e-v[i]) <= eps";
  output Integer result[20]
    "v[result] = e (first occurrence of e); result=0, if not found";

protected
  Integer i;
  Integer k;
algorithm
  k := 0;
  i := 1;
  while i <= size(v, 1) loop
    if abs(v[i] - e) <= eps then
      k :=k + 1;
      result[k] := i;
      i := i + 1;
    else
      i := i + 1;
    end if;
  end while;

   assert(k<size(result,1),"Need to increase maximum size of result vector");
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Vectors.<b>find</b>(e, v);
Vectors.<b>find</b>(e, v, eps=0);
</pre></blockquote>
<h4>Description</h4>
<p>
The function call \"<code>Vectors.find(e, v)</code>\" returns the index of the first occurrence of input e in vector <b>v</b>.
The test of equality is performed by \"abs(e-v[i]) &le; eps\", where \"eps\"
can be provided as third argument of the function. Default is \"eps = 0\".
</p>
<h4>Example</h4>
<blockquote><pre>
  Real v[3] = {1, 2, 3};
  Real e1 = 2;
  Real e2 = 3.01;
  Boolean result;
<b>algorithm</b>
  result := Vectors.find(e1,v);          // = <b>2</b>
  result := Vectors.find(e2,v);          // = <b>0</b>
  result := Vectors.find(e2,v,eps=0.1);  // = <b>3</b>
</pre></blockquote>
<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Vectors.isEqual\">Vectors.isEqual</a>
</p>
</html>"));
end find_all;
