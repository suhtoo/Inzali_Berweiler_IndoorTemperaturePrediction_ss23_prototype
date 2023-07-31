package BaseLib
  class HeatCond
    extends BaseLib.TwoPort;
    parameter Modelica.SIunits.Area A = 16;
    parameter Modelica.SIunits.Thickness d = 0.1;
    parameter Modelica.SIunits.ThermalConductivity lambda = 2.4;
  equation
    // no storage of heat
    port_a.Q_flow + port_b.Q_flow = 0;
    port_a.Q_flow = lambda * A / d * (port_a.T - port_b.T);
    annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Rectangle(extent = {{60, 60}, {80, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {244, 244, 244}, fillPattern = FillPattern.Solid), Rectangle(extent = {{40, 60}, {60, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {216, 216, 216}, fillPattern = FillPattern.Solid), Rectangle(extent = {{20, 60}, {40, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {207, 207, 207}, fillPattern = FillPattern.Solid), Rectangle(extent = {{0, 60}, {20, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {188, 188, 188}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-20, 60}, {0, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {182, 182, 182}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-40, 60}, {-20, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {173, 173, 173}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-60, 60}, {-40, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {168, 168, 168}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-80, 60}, {-60, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {156, 156, 156}, fillPattern = FillPattern.Solid)}), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Rectangle(extent = {{60, 60}, {80, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {244, 244, 244},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Rectangle(extent = {{40, 60}, {60, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {216, 216, 216},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Rectangle(extent = {{20, 60}, {40, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {207, 207, 207},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Rectangle(extent = {{0, 60}, {20, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {188, 188, 188},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Rectangle(extent = {{-20, 60}, {0, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {182, 182, 182},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Rectangle(extent = {{-40, 60}, {-20, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {173, 173, 173},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Rectangle(extent = {{-60, 60}, {-40, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {168, 168, 168},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Rectangle(extent = {{-80, 60}, {-60, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {156, 156, 156},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid)}), Window(x = 0.4, y = 0.4, width = 0.6, height = 0.6), Documentation(revisions = "<html>
<p><ul>
<li><i>April 10, 2013&nbsp;</i> by Ole Odendahl:<br/>Added documentation and formatted appropriately</li>
</ul></p>
</html>", info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>The <b>HeatCond</b> model represents the phenomenon of heat conduction.</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
</html>"));
  end HeatCond;

  class Load "Heat load"
    extends BaseLib.OnePort;
    parameter Modelica.SIunits.Density rho = 1600;
    parameter Modelica.SIunits.SpecificHeatCapacity c = 1000;
    parameter Modelica.SIunits.Thickness d = 0.2;
    parameter Modelica.SIunits.Area A = 16;
    Modelica.SIunits.Mass m;
  equation
    m = rho * A * d;
    der(T) = 1 / m / c * port_a.Q_flow;
    annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}, lineThickness = 0.5, fillPattern = FillPattern.Sphere, fillColor = {191, 191, 191})}), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0},
              lineThickness =                                                                                                                                                                                                        0.5,
              fillPattern =                                                                                                                                                                                                        FillPattern.Sphere, fillColor = {227, 227, 227})}), Window(x = 0.3, y = 0.18, width = 0.6, height = 0.6), Documentation(revisions = "<html>
<p><ul>
<li><i>April 10, 2013&nbsp;</i> by Ole Odendahl:<br/>Added documentation and formatted appropriately</li>
</ul></p>
</html>", info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>The <b>Load</b> model represents a heat capacity, which is described by its area, density, thickness and material specific heat capacity.</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
</html>"));
  end Load;

  partial class OnePort
    parameter Modelica.SIunits.Temperature T0 = Modelica.SIunits.Conversions.from_degC(16) "Initial temperature";
    Modelica.SIunits.Temperature T(start = T0) "Initial temperature";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(Placement(visible = true, transformation(origin = {-97, 3}, extent = {{-15, -15}, {15, 15}}, rotation = 0), iconTransformation(origin = {-100, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    T = port_a.T;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end OnePort;

  partial class TwoPort
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(Placement(visible = true, transformation(origin = {-97, 1}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {-100, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation(Placement(visible = true, transformation(origin = {100, 8.88178e-16}, extent = {{-22, -22}, {22, 22}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end TwoPort;

  class I_to_Qdot "Compute the heat flow caused by radiation on a surface"
    parameter Modelica.SIunits.Area A = 2 "Area of surface";
    parameter Modelica.SIunits.Area frameratio = 0 "frame ratio if window";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(Placement(visible = true, transformation(origin = {101, 1}, extent = {{-23, -23}, {23, 23}}, rotation = 0), iconTransformation(origin = {90, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.RadiationPort radiationport1 annotation(Placement(visible = true, transformation(origin = {-98, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-93, 7}, extent = {{-15, -21}, {15, 9}}, rotation = 0)));
  equation
    port_a.Q_flow = -radiationport1.I * A * (1 - frameratio);
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics={  Rectangle(origin = {0.17, -10.73}, fillColor = {103, 103, 103},
              fillPattern =                                                                                                                                                                                                        FillPattern.HorizontalCylinder, extent = {{-70.13, 76.73}, {69.8, -59.57}}), Text(origin = {0.04, 2.02}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{-59.9, 46.37}, {59.9, -46.37}}, textString = "I to Qdot")}));
  end I_to_Qdot;

  connector Input_Radiation "Scalar total radiation connector (output)"
    input Modelica.SIunits.RadiantEnergyFluenceRate I;
    annotation(Documentation(info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>
The <b>oc_total_rad</b> connector is used for scalar total radiation output.
</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
</html>
", revisions = "<html>
<p><ul>
<li><i>April 10, 2013&nbsp;</i> by Ole Odendahl:<br/>Formatted documentation appropriately</li>
</ul></p>
</html>"), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics={  Rectangle(origin = {159.608, 0.190678}, fillColor = {215, 215, 215},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{-100, 100}, {-60, -100}}), Polygon(origin = {-0.792079, 0.858086}, fillColor = {255, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, points = {{-60, 6}, {-60, 28}, {-60, 28}, {-8, 0}, {-60, -26}, {-60, -26}, {-60, -6}, {-60, 6}}), Line(origin = {-0.792079, 0.858086}, points = {{-40, -80}, {-40, -62}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{20, -22}, {74, -40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{8, -44}, {46, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-14, -58}, {-2, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-14, 58}, {-2, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{20, 22}, {72, 40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{6, 44}, {46, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{24, 0}, {78, 0}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-40, 62}, {-40, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Ellipse(origin = {-0.792079, 0.858086}, lineColor = {255, 128, 0}, extent = {{18, 58}, {-98, -58}}, endAngle = 360)}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics = {Ellipse(lineColor = {255, 128, 0}, extent = {{18, 58}, {-98, -58}}, endAngle = 360), Line(points = {{-40, 62}, {-40, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{24, 0}, {78, 0}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{6, 44}, {46, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{20, 22}, {72, 40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-14, 58}, {-2, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-14, -58}, {-2, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{8, -44}, {46, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{20, -22}, {74, -40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-40, -80}, {-40, -62}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Rectangle(origin = {159.862, 0}, fillColor = {215, 215, 215}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {-60, -100}}), Polygon(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, points = {{-60, 6}, {-60, 28}, {-60, 28}, {-8, 0}, {-60, -26}, {-60, -26}, {-60, -6}, {-60, 6}})}));
  end Input_Radiation;

  model TwoStar_RadEx "Adaptor for approximative longwave radiation exchange"
    parameter Modelica.SIunits.Area A = 12 "Area of radiation exchange";
    parameter Modelica.SIunits.Emissivity eps = 0.95 "Emissivity";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(Placement(visible = true, transformation(origin = {-96, -1}, extent = {{-14, -15}, {14, 15}}, rotation = 0), iconTransformation(origin = {-100, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.Star star annotation(Placement(visible = true, transformation(origin = {98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    port_a.Q_flow + star.Q_flow = 0;
    port_a.Q_flow = Modelica.Constants.sigma * eps * A * (port_a.T * port_a.T * port_a.T * port_a.T - star.T * star.T * star.T * star.T);
    annotation(Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {135, 150, 177}, fillPattern = FillPattern.Solid), Text(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 0, 0}, pattern = LinePattern.None, fillColor = {135, 150, 177}, fillPattern = FillPattern.Solid, textString = "2*")}), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {135, 150, 177},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Text(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 0, 0}, pattern = LinePattern.None, fillColor = {135, 150, 177},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, textString = "2*")}), Documentation(info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>The <b>TwoStar_RadEx</b> model cobines the <b><a href=\"Interfaces.port_a\">port_a</a></b> connector and the <b><a href=\"Interfaces.Star\">Star</a></b> connector. To model longwave radiation exchange of a surfaces, just connect the <b><a href=\"Interfaces.port_a\">port_a</a></b> connector to the outmost layer of the surface and connect the <b><a href=\"Interfaces.Star\">Star</a></b> connector to the <b><a href=\"Interfaces.Star\">Star</a></b> connectors of an unlimited number of corresponding surfaces. </p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
<p><h4><font color=\"#008000\">Concept</font></h4></p>
<p>Since exact calculation of longwave radiation exchange inside a room demands for the computation of view factors, it may be very complex to achieve for non-rectangular room layouts. Therefore, an approximated calculation of radiation exchange basing on the proportions of the affected surfaces is an alternative. The underlying concept of this approach is known as the &QUOT;two star&QUOT; room model. </p>
</html>", revisions = "<html>
<ul>
  <li><i>April 10, 2013&nbsp;</i> by Ole Odendahl:<br/>Formatted documentation appropriately</li>
  <li><i>June 16, 2006&nbsp;</i>
         by Timo Haase:<br>
         Implemented.</li>
</ul>
</html>"));
  end TwoStar_RadEx;

  connector Star "Connector for twostar (approximated) radiation exchange"
    Modelica.SIunits.Temperature T;
    flow Modelica.SIunits.HeatFlowRate Q_flow;
    annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-80, 80}, {80, -80}}, lineColor = {95, 95, 95}, pattern = LinePattern.None, fillColor = {255, 255, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Polygon(points = {{-13, 86}, {13, 86}, {13, 12}, {77, 34}, {85, 6}, {22, -14}, {62, -72}, {37, -88}, {0, -28}, {-35, -88}, {-60, -72}, {-22, -14}, {-85, 6}, {-77, 34}, {-13, 12}, {-13, 86}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-82, 84}, {78, -76}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {255, 255, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Polygon(points = {{-13, 86}, {13, 86}, {13, 12}, {77, 34}, {85, 6}, {22, -14}, {62, -72}, {37, -88}, {0, -28}, {-35, -88}, {-60, -72}, {-22, -14}, {-85, 6}, {-77, 34}, {-13, 12}, {-13, 86}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid)}), Documentation(info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>
The <b>Star</b> connector is technically a clone of the <a href=Interfaces.Therm><b>Therm</b></a> connector. But the carried data has to be interpreted in a different way: the temperature T is a virtual temperature describing the potential of longwave radiation exchange inside the room. The heat flow Q_flow is the resulting energy flow due to longwave radiation.
</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
</html>", revisions = "<html>
<ul>
<li><i>April 10, 2013&nbsp;</i> by Ole Odendahl:<br/>Formatted documentation appropriately</li>
<li><i>July 12, 2009&nbsp;</i>
         by Peter Matthes:<br>
         Switched to Modelica.SIunits.Temperature.</li>
  <li><i>June 16, 2006&nbsp;</i>
         by Timo Haase:<br>
         Implemented.</li>
  
</ul>
</html>"));
  end Star;

  model Waermeleiter "Lumped thermal element transporting heat without storing it"
    extends Modelica.Thermal.HeatTransfer.Interfaces.Element1D;
    parameter Modelica.SIunits.CoefficientOfHeatTransfer U "Constant thermal conductivity of material";
    parameter Modelica.SIunits.Area A "area";
  equation
    Q_flow = U * A * dT;
    annotation(Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-90, 70}, {90, -70}}, lineColor = {0, 0, 0}, pattern = LinePattern.None, fillColor = {192, 192, 192},
              fillPattern =                                                                                                                                                                                                        FillPattern.Backward), Line(points = {{-90, 70}, {-90, -70}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{90, 70}, {90, -70}}, color = {0, 0, 0}, thickness = 0.5), Text(extent = {{-150, 115}, {150, 75}}, textString = "%name", lineColor = {0, 0, 255}), Text(extent = {{-150, -75}, {150, -105}}, lineColor = {0, 0, 0}, textString = "U=%U,A=%A")}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Line(points = {{-80, 0}, {80, 0}}, color = {255, 0, 0}, thickness = 0.5, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{-100, -20}, {100, -40}}, lineColor = {255, 0, 0}, textString = "Q_flow"), Text(extent = {{-100, 40}, {100, 20}}, lineColor = {0, 0, 0}, textString = "dT = port_a.T - port_b.T")}), Documentation(info = "<HTML>
<p>
This is a model for transport of heat without storing it; see also:
<a href=\"modelica://Modelica.Thermal.HeatTransfer.Components.ThermalResistor\">ThermalResistor</a>.
It may be used for complicated geometries where
the thermal conductance G (= inverse of thermal resistance)
is determined by measurements and is assumed to be constant
over the range of operations. If the component consists mainly of
one type of material and a regular geometry, it may be calculated,
e.g., with one of the following equations:
</p>
<ul>
<li><p>
    Conductance for a <b>box</b> geometry under the assumption
    that heat flows along the box length:</p>
    <pre>
    G = k*A/L
    k: Thermal conductivity (material constant)
    A: Area of box
    L: Length of box
    </pre>
    </li>
<li><p>
    Conductance for a <b>cylindrical</b> geometry under the assumption
    that heat flows from the inside to the outside radius
    of the cylinder:</p>
    <pre>
    G = 2*pi*k*L/log(r_out/r_in)
    pi   : Modelica.Constants.pi
    k    : Thermal conductivity (material constant)
    L    : Length of cylinder
    log  : Modelica.Math.log;
    r_out: Outer radius of cylinder
    r_in : Inner radius of cylinder
    </pre>
    </li>
</ul>
<pre>
    Typical values for k at 20 degC in W/(m.K):
      aluminium   220
      concrete      1
      copper      384
      iron         74
      silver      407
      steel        45 .. 15 (V2A)
      wood         0.1 ... 0.2
</pre>
</html>"));
  end Waermeleiter;

  connector TwoRadiationPort "radiation connector"
    input Modelica.SIunits.RadiantEnergyFluenceRate I_dir;
    input Modelica.SIunits.RadiantEnergyFluenceRate I_diff;
    annotation(Documentation(info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>
The <b>oc_total_rad</b> connector is used for scalar total radiation output.
</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
</html>
", revisions = "<html>
<p><ul>
<li><i>April 10, 2013&nbsp;</i> by Ole Odendahl:<br/>Formatted documentation appropriately</li>
</ul></p>
</html>"), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics={  Rectangle(origin = {159.608, 0.190678}, fillColor = {215, 215, 215},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{-100, 100}, {-60, -100}}), Polygon(origin = {-2.79208, 22.8581}, fillColor = {255, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, points = {{-60, 6}, {-60, 28}, {-60, 28}, {-8, 0}, {-60, -26}, {-60, -26}, {-60, -6}, {-60, 6}}), Line(origin = {-0.792079, 0.858086}, points = {{-40, -80}, {-40, -62}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{20, -22}, {74, -40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{8, -44}, {46, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-14, -58}, {-2, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-14, 58}, {-2, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{20, 22}, {72, 40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{6, 44}, {46, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{24, 0}, {78, 0}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-40, 62}, {-40, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Ellipse(origin = {-0.792079, 0.858086}, lineColor = {255, 128, 0}, extent = {{18, 58}, {-98, -58}}, endAngle = 360), Polygon(origin = {-2.79208, -17.1419}, fillColor = {255, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, points = {{-60, 6}, {-60, 28}, {-60, 28}, {-8, 0}, {-60, -26}, {-60, -26}, {-60, -6}, {-60, 6}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics={  Ellipse(lineColor = {255, 128, 0}, extent = {{18, 58}, {-98, -58}}, endAngle = 360), Line(points = {{-40, 62}, {-40, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{24, 0}, {78, 0}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{6, 44}, {46, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{20, 22}, {72, 40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-14, 58}, {-2, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-14, -58}, {-2, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{8, -44}, {46, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{20, -22}, {74, -40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-40, -80}, {-40, -62}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Rectangle(origin = {159.862, 0}, fillColor = {215, 215, 215},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{-100, 100}, {-60, -100}}), Polygon(origin = {0, 18}, fillColor = {255, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, points = {{-60, 6}, {-60, 28}, {-60, 28}, {-8, 0}, {-60, -26}, {-60, -26}, {-60, -6}, {-60, 6}}), Polygon(origin = {0, -16}, fillColor = {255, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, points = {{-60, 6}, {-60, 28}, {-60, 28}, {-8, 0}, {-60, -26}, {-60, -26}, {-60, -6}, {-60, 6}})}));
  end TwoRadiationPort;

  connector RadiationPort "radiation connector"
    input Modelica.SIunits.RadiantEnergyFluenceRate I;
    annotation(Documentation(info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>
The <b>oc_total_rad</b> connector is used for scalar total radiation output.
</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
</html>
", revisions = "<html>
<p><ul>
<li><i>April 10, 2013&nbsp;</i> by Ole Odendahl:<br/>Formatted documentation appropriately</li>
</ul></p>
</html>"), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics={  Rectangle(origin = {159.608, 0.190678}, fillColor = {215, 215, 215},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{-100, 100}, {-60, -100}}), Polygon(origin = {-2.79208, 0.858101}, fillColor = {255, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, points = {{-60, 6}, {-60, 28}, {-60, 28}, {-8, 0}, {-60, -26}, {-60, -26}, {-60, -6}, {-60, 6}}), Line(origin = {-0.792079, 0.858086}, points = {{-40, -80}, {-40, -62}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{20, -22}, {74, -40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{8, -44}, {46, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-14, -58}, {-2, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-14, 58}, {-2, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{20, 22}, {72, 40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{6, 44}, {46, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{24, 0}, {78, 0}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(origin = {-0.792079, 0.858086}, points = {{-40, 62}, {-40, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Ellipse(origin = {-0.792079, 0.858086}, lineColor = {255, 128, 0}, extent = {{18, 58}, {-98, -58}}, endAngle = 360)}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics={  Ellipse(lineColor = {255, 128, 0}, extent = {{18, 58}, {-98, -58}}, endAngle = 360), Line(points = {{-40, 62}, {-40, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{24, 0}, {78, 0}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{6, 44}, {46, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{20, 22}, {72, 40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-14, 58}, {-2, 80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-14, -58}, {-2, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{8, -44}, {46, -80}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{20, -22}, {74, -40}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Line(points = {{-40, -80}, {-40, -62}}, color = {255, 128, 0}, smooth = Smooth.Bezier), Rectangle(origin = {159.862, 0}, fillColor = {215, 215, 215},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{-100, 100}, {-60, -100}}), Polygon(fillColor = {255, 0, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, points = {{-60, 6}, {-60, 28}, {-60, 28}, {-8, 0}, {-60, -26}, {-60, -26}, {-60, -6}, {-60, 6}})}));
  end RadiationPort;

  connector TwoRealPort = input Real "'input Real' as connector" annotation(defaultComponentName = "u", Icon(graphics={  Polygon(lineColor = {0, 0, 127}, fillColor = {0, 0, 127},
            fillPattern =                                                                                                                                                                        FillPattern.Solid, points = {{-100.0, 100.0}, {100.0, 0.0}, {-100.0, -100.0}})}, coordinateSystem(extent = {{-100.0, -100.0}, {100.0, 100.0}}, preserveAspectRatio = true, initialScale = 0.2)), Diagram(coordinateSystem(preserveAspectRatio = true, initialScale = 0.2, extent = {{-100.0, -100.0}, {100.0, 100.0}}), graphics = {Polygon(lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid, points = {{0.0, 50.0}, {100.0, 0.0}, {0.0, -50.0}, {0.0, 50.0}}), Text(lineColor = {0, 0, 127}, extent = {{-10.0, 60.0}, {-10.0, 85.0}}, textString = "%name")}), Documentation(info = "<html>
                           <p>
                           Connector with one input signal of type Real.
                           </p>
                           </html>"));

  connector SolarPositionPort "connector for declination and hour angle"
    input Real omega "hour angle";
    input Real delta "declination";
    annotation(defaultComponentName = "solarposition", Documentation(info = "<html>
 <p>
 Connector with one input signal of type Real.
 </p>
 </html>"), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2, grid = {2, 2}), graphics={  Polygon(lineColor = {0, 0, 127}, fillColor = {255, 255, 0},
              fillPattern =                                                                                                                                                                                                     FillPattern.Solid, points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}), Text(lineColor = {0, 0, 127}, extent = {{-10, 60}, {-10, 85}}, textString = "%name")}), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.2, grid = {2, 2}), graphics={  Polygon(lineColor = {0, 0, 127}, fillColor = {255, 255, 0},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}})}));
  end SolarPositionPort;

  model addRadiation
    parameter Real alpha_SW = 0.7 "shortwave absorption outside";
    BaseLib.TwoRadiationPort tworadiationport1 annotation(Placement(visible = true, transformation(origin = {-98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.RadiationPort radiationport1 annotation(Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    radiationport1.I = alpha_SW*(tworadiationport1.I_diff + tworadiationport1.I_dir);
    annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics={  Line(origin = {0, 17.5}, points = {{0, 30.5}, {0, -51.5}, {0, -49.5}}, thickness = 20), Line(origin = {-44, 49.5}, points = {{102, -51.5}, {0, -51.5}, {-10, -51.5}}, thickness = 20)}));
  end addRadiation;

  class Airload "Air volume"
    extends BaseLib.OnePort;
    parameter Modelica.SIunits.Density rho = 1.19 "Density of air";
    parameter Modelica.SIunits.SpecificHeatCapacity c = 1007 "Specific heat capacity of air";
    parameter Modelica.SIunits.Volume V = 48.0 "Volume of the room";
    Modelica.SIunits.Mass m;
    Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_degC;
  equation
    m = rho * V;
    der(T) = 1 / m / c * port_a.Q_flow;
    T_degC = Modelica.SIunits.Conversions.to_degC(T);
    annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}, fillColor = {211, 243, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-28, 14}, {32, -52}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, textString = "Air")}), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics={  Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}, fillColor = {211, 243, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Text(extent = {{-30, 16}, {30, -50}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, textString = "Air")}), Window(x = 0.25, y = 0.09, width = 0.6, height = 0.6), Documentation(revisions = "<html>
<p><ul>
<li><i>May 02, 2013&nbsp;</i> by Ole Odendahl:<br/>Formatted documentation appropriately</li>
</ul></p>
</html>", info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>
The <b>Airload</b> model represents a heat capacity consisting of air. It is described by its volume, density and specific heat capacity.
</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
</html>"));
  end Airload;

  model Integration
    TwoRadiationPort in_I annotation(Placement(visible = true, transformation(origin = {-52, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  end Integration;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
end BaseLib;
