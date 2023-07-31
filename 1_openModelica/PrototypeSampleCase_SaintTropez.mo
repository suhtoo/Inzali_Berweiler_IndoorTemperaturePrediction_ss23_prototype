package HomeAssignment_2023
  model window "WindowModel with Ventilation"
    parameter Real A = 2.0 "window area";
    parameter Real frameratio = 0.2 "window frame ration";
    parameter Real N_pane = 1.526 "refraction index/Brechungsindex";
    parameter Real d = 0.004 "pane thickness";
    parameter Integer number_panes = 2 "number of panes in the window";
    parameter Real beta = 90 "tilt angle";
    parameter Real gamma = 0 "azimuth angle";
    parameter Real rho_floor = 0.2 "ground reflection coefficient";
    parameter Real K = 4.0 "extinction parameter";
    parameter Real U = 2.1 "U value";
    parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg eps = 50.789389 "latitude of location";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatport_outside annotation(
      Placement(visible = true, transformation(origin = {-100, -72}, extent = {{-14, -14}, {14, 14}}, rotation = 0), iconTransformation(origin = {-102, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatport_inside annotation(
      Placement(visible = true, transformation(origin = {98, -72}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {100, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.TwoRadiationPort tworadiationport1 annotation(
      Placement(visible = true, transformation(origin = {-96, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-116, 8}, {-96, 28}}, rotation = 0)));
    BaseLib.TwoStar_RadEx twostar_radex1(A = 1, eps = 0.95) annotation(
      Placement(visible = true, transformation(origin = {70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.Star star1 annotation(
      Placement(visible = true, transformation(origin = {102, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.SolarPositionPort solarposition annotation(
      Placement(visible = true, transformation(origin = {-100, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-94, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.GeneigteFlaeche tiltedSurface1(beta = beta, eps = eps, gamma = gamma, rho_floor = rho_floor) annotation(
      Placement(visible = true, transformation(origin = {-53, 59}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
    BaseLib.HeatCond heatConductor1(A = A, d = d, lambda = U * d) annotation(
      Placement(visible = true, transformation(origin = {-7, -71}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    HomeAssignment_2023.RadiativeHeatTransfer radiative_heat_transfer1(N_pane = N_pane, K = K, d = d, number_panes = number_panes) annotation(
      Placement(visible = true, transformation(origin = {2, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    HomeAssignment_2023.ventilation ventilation1(V = 76.8) annotation(
      Placement(visible = true, transformation(origin = {-5, -11}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    Modelica.Blocks.Interfaces.IntegerInput u annotation(
      Placement(visible = true, transformation(origin = {-100, -18}, extent = {{-14, -14}, {14, 14}}, rotation = 0), iconTransformation(origin = {-98, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.I_to_Qdot i_to_Qdot(A = A, frameratio = frameratio) annotation(
      Placement(transformation(extent = {{32, 50}, {52, 70}})));
  equation
    connect(ventilation1.WindowStatus, u) annotation(
      Line(points = {{-23.62, -2.64}, {-78, -2.64}, {-78, -18}, {-100, -18}, {-100, -18}}, color = {255, 127, 0}));
    connect(ventilation1.port_out, heatport_outside) annotation(
      Line(points = {{-24, -19.36}, {-50, -19.36}, {-50, -72}, {-100, -72}, {-100, -72}}, color = {191, 0, 0}));
    connect(ventilation1.port_in, heatport_inside) annotation(
      Line(points = {{14, -19.36}, {50, -19.36}, {50, -70}, {98, -70}, {98, -72}}, color = {191, 0, 0}));
    connect(tiltedSurface1.Out_phi, radiative_heat_transfer1.In_phi_1) annotation(
      Line(points = {{-34.1, 43.88}, {-30, 43.88}, {-30, 50}, {-18, 50}, {-18, 49.6}}, color = {0, 0, 127}));
    connect(tiltedSurface1.out_rad, radiative_heat_transfer1.tworadiationport1) annotation(
      Line(points = {{-35.36, 63.62}, {-32, 63.62}, {-32, 74}, {-18, 74}, {-18, 73.6}}));
    connect(heatConductor1.port_b, heatport_inside) annotation(
      Line(points = {{10, -71}, {98, -71}, {98, -72}, {98, -72}}, color = {191, 0, 0}));
    connect(heatport_outside, heatConductor1.port_a) annotation(
      Line(points = {{-100, -72}, {-24, -72}, {-24, -70.66}, {-24, -70.66}}, color = {191, 0, 0}));
    connect(solarposition, tiltedSurface1.solarposition) annotation(
      Line(points = {{-100, 66}, {-72, 66}, {-72, 65.72}, {-71.9, 65.72}}, color = {0, 0, 127}));
    connect(tworadiationport1, tiltedSurface1.in_rad) annotation(
      Line(points = {{-96, 26}, {-82, 26}, {-82, 44}, {-73.16, 44}, {-73.16, 44.3}}, color = {255, 128, 0}));
    connect(twostar_radex1.star, star1) annotation(
      Line(points = {{79.8, 60}, {102, 60}}, color = {95, 95, 95}));
    connect(radiative_heat_transfer1.radiationport1, i_to_Qdot.radiationport1) annotation(
      Line(points = {{21.6, 60}, {27.15, 60}, {27.15, 60.1}, {32.7, 60.1}}, color = {0, 0, 0}));
    connect(i_to_Qdot.port_a, twostar_radex1.port_a) annotation(
      Line(points = {{51, 59.8}, {58, 59.8}, {58, 60.4}, {60, 60.4}}, color = {191, 0, 0}));
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {2.88765, -3.00531}, lineColor = {170, 85, 0}, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 3, extent = {{-80.76000000000001, 86.65000000000001}, {80.76000000000001, -86.65000000000001}}), Rectangle(origin = {-95.2436, 64.82559999999999}, lineColor = {170, 85, 0}, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 3, extent = {{16.6816, 18.7968}, {100.782, -154.503}}), Rectangle(origin = {-95.2659, 64.8034}, lineColor = {170, 85, 0}, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 3, extent = {{16.6816, 18.7968}, {179.091, -32.3675}})}),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end window;

  model ventilation
    parameter Real rho = 1.204 "air density in kg/m^3";
    parameter Real n_open = 20 "air change rate for open window";
    parameter Real n_tilt = 2 "air change rate for tilted window";
    parameter Real n_50 = 3 "air change rate at 50 pa pressure difference";
    parameter Real e = 0.05 "wind shielding factor";
    parameter Real epsilon = 1 "height correction factor";
    parameter Real V = 60 "room volume in m^3";
    parameter Real c_p = 1005 "specific heat capicity of air in J/K*kg";
    Real delta_T "temperature difference in K";
    Real V_dot "volume flow in m^3/s";
    Real Q_flow "heat flow in W";
    Modelica.Blocks.Interfaces.IntegerInput WindowStatus annotation(
      Placement(visible = true, transformation(origin = {-104, 44}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-98, 44}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_out annotation(
      Placement(visible = true, transformation(origin = {-100, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_in annotation(
      Placement(visible = true, transformation(origin = {100, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    if WindowStatus == 2 then
      V_dot = n_tilt * V / 3600;
    elseif WindowStatus == 1 then
      V_dot = n_50 * e * epsilon * 2 * V / 3600;
    else
      V_dot = n_open * V / 3600;
    end if;
    port_in.Q_flow + port_out.Q_flow = 0;
    Q_flow = port_in.Q_flow;
    delta_T = port_in.T - port_out.T;
    Q_flow = V_dot * c_p * rho * delta_T;
    annotation(
      Icon(graphics = {Rectangle(origin = {-4, 4}, fillColor = {170, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-74, 80}, {74, -80}}), Line(origin = {9.92, -3.92}, points = {{-43.916, 1.9238}, {30.084, 1.9238}, {30.084, 15.9238}, {44.084, 1.9238}, {32.084, -16.0762}, {32.084, -6.0762}, {-43.916, -6.0762}, {-43.916, 1.9238}, {-43.916, 1.9238}}, thickness = 3.5), Line(origin = {-0.08, 28.08}, rotation = 180, points = {{-43.916, 1.9238}, {30.084, 1.9238}, {30.084, 15.9238}, {44.084, 1.9238}, {32.084, -16.0762}, {32.084, -6.0762}, {-43.916, -6.0762}, {-43.916, 1.9238}, {-43.916, 1.9238}}, thickness = 3.5)}));
  end ventilation;

  model Input_Window_Status
    Modelica.Blocks.Interfaces.IntegerOutput Window_Status annotation(
      Placement(visible = true, transformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    if time <= 1 then
      Window_Status = 2;
    else
      Window_Status = 1;
    end if;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Input_Window_Status;

  model N_Layer
    extends BaseLib.TwoPort;
    parameter Modelica.SIunits.Height h = 3 "Height" annotation(
      Dialog(group = "Geometry"));
    parameter Modelica.SIunits.Length l = 4 "Length" annotation(
      Dialog(group = "Geometry"));
    parameter Modelica.SIunits.Area clearance = 0 "Area of clearance" annotation(
      Dialog(group = "Geometry"));
    parameter Integer n(min = 1) = 8 "Number of layers" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    // to be re-implemented
    parameter Modelica.SIunits.Thickness d[n] = fill(0.1, n) "Thickness" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    // parameter  Modelica.SIunits.Thickness d [:] = {0.1, 0.1, 0.1} "Thickness" annotation(
    //   Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.Density rho[n] = fill(1600, n) "Density" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.ThermalConductivity lambda[n] = fill(2.4, n) "Thermal conductivity" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.SpecificHeatCapacity c[n] = fill(1000, n) "Specific heat capacity" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.Temperature T0 = Modelica.SIunits.Conversions.from_degC(16) "Initial temperature" annotation(
      Dialog(group = "Thermal"));
    // 2n HeatConds
    BaseLib.HeatCond HeatConda[n](each A = A, d = d / 2, lambda = lambda, port_a(each T(start = T0)), port_b(each T(start = T0))) annotation(
      Placement(transformation(extent = {{8, -8}, {28, 12}}, rotation = 0)));
    BaseLib.HeatCond HeatCondb[n](each A = A, d = d / 2, lambda = lambda, port_a(each T(start = T0)), port_b(each T(start = T0))) annotation(
      Placement(transformation(extent = {{-50, -8}, {-30, 12}}, rotation = 0)));
    // n Loads
    BaseLib.Load Load[n](each T0 = T0, each A = A, d = d, rho = rho, c = c) annotation(
      Placement(transformation(extent = {{-8, -62}, {12, -42}}, rotation = 0)));
  protected
    parameter Modelica.SIunits.Area A = h * l - clearance;
  equation
// connecting inner elements HeatCondb[i]--Load[i]--HeatConda[i] to n groups
    connect(HeatCondb[1].port_b, Load[1].port_a) annotation(
      Line(points = {{-30, 2}, {-20, 2}, {-20, -52}, {-8, -52}, {-8, -51.8}}, color = {191, 0, 0}));
    connect(HeatConda[1].port_a, Load[1].port_a) annotation(
      Line(points = {{8, 2.2}, {-20, 2.2}, {-20, -52}, {-8, -52}, {-8, -51.8}, {-8, -51.8}}, color = {191, 0, 0}));
    for i in 2:n loop
      connect(HeatCondb[i].port_b, Load[i].port_a);
      connect(HeatConda[i].port_a, Load[i].port_a);
    end for;
// establishing n-1 connections of HeatCondb--Load--HeatConda groups
    for i in 1:n - 1 loop
      connect(HeatConda[i].port_b, HeatCondb[i + 1].port_a);
    end for;
// connecting outmost elements to connectors: port_a--HeatCondb[1]...HeatConda[n]--HeatConv1--port_b
    connect(HeatCondb[1].port_a, port_a) annotation(
      Line(points = {{-50, 2.2}, {-69.5, 2.2}, {-69.5, 1}, {-97, 1}}, color = {200, 100, 0}));
    connect(HeatConda[n].port_b, port_b) annotation(
      Line(points = {{28, 2}, {98, 2}, {98, 0}, {100, 0}}, color = {191, 0, 0}));
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {-2.79, 0}, fillColor = {170, 62, 54}, fillPattern = FillPattern.Vertical, lineThickness = 2, extent = {{-43.07, 98.02}, {43.07, -98.02}}), Text(origin = {-35.86, 75.67}, extent = {{0, -0.23}, {64.9593, -44.2403}}, textString = "Wall"), Line(origin = {-69.03, -0.7}, points = {{-21.3038, 1.11022e-016}, {21.3038, 1.11022e-016}}, thickness = 1), Line(origin = {64.81999999999999, 0.36}, points = {{-22.6808, -1.06181}, {25.5148, -1.06181}, {17.3658, -0.828977}}, thickness = 1)}),
      Documentation(info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>The <b>ConvNLayerClearanceStar</b> model represents a wall, consisting of n different layers with natural convection on one side and (window) clearance.</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
<p><h4><font color=\"#008000\">Concept</font></h4></p>
<p>There is one inner and one outer <b><a href=\"BaseLib.Therm\">Therm</a></b>-connector to simulate one-dimensional heat transfer through the wall and heat storage within the wall.</p>
<p>The <b>ConvNLayerClearanceStar</b> model extends the basic concept by adding the functionality of approximated longwave radiation exchange. Simply connect all radiation exchanging surfaces via their <b><a href=\"Modelica://BaseLib.Star\">Star</a></b>-connectors. </p>
<p><b><font style=\"color: #ff0000; \">Attention:</font></b> The first element in each vector represents the layer connected to <code>port_a</code>, the last element represents the layer connected to <code>port_b</code>. </p>
<p><h4><font color=\"#008000\">Example Results</font></h4></p>
<p>This model is part of <a href=\"Building.Components.Walls.WallSimple\">WallSimple</a> or <a href=\"Building.Components.Walls.WallSimple_newconnectors\">WallSimple_newconnectors</a> and therefore also part of the corresponding examples <a href=\"Building.Examples.Walls.InsideWall\">InsideWall</a> and <a href=\"Building.Examples.Walls.OutsideWall\">OutsideWall</a>. </p>
</html>", revisions = "<html>
<ul>
<li><i>May 02, 2013&nbsp;</i> by Ole Odendahl:<br/>Formatted documentation appropriately</li>
<li><i>Aug. 08, 2006&nbsp;</i>
       by Peter Matthes:<br>
       Fixed wrong connection with heatConv-Module and added connection graphics.</li>

<li><i>June 19, 2006&nbsp;</i>
       by Timo Haase:<br>
       Implemented.</li>
</ul>
</html>"),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {0, 12}, extent = {{-80, 60}, {80, -100}})}));
  end N_Layer;

  class HeatConvUniversal
    extends BaseLib.TwoPort;
    constant Integer hor_up = 1 "Horizontal surface top side";
    constant Integer hor_down = 2 "Horizontal surface down side";
    constant Integer vert_interior = 3 "Heat convection at vertical interior wall";
    constant Integer vert_exterior = 4 "Heat convection at outside wall";
    constant Integer custom = 5 "Custom set coefficient of heat transfer";

    type Temp
      extends Integer;
      annotation(
        Evaluate = true,
        choices(choice = hor_up "Horizontal surface top side", choice = hor_down "Horizontal surface down side", choice = vert_interior "Heat convection at vertical interior wall", choice = vert_exterior "Heat convection at outside wall", choice = custom "Custom set coefficient of heat transfer"));
    end Temp;

    parameter Temp control_type = vert_exterior annotation(
      Dialog(group = "Control"));
    Modelica.SIunits.CoefficientOfHeatTransfer alpha;
    parameter Modelica.SIunits.Area A = 16 "Area of surface";
    parameter Modelica.SIunits.Velocity v = 3 "Outside air velocity [m/s], if not given via WindSpeedPort" annotation(
      Dialog(enable = control_type == vert_exterior));
    parameter Modelica.SIunits.CoefficientOfHeatTransfer alpha_custom = 2 "Coefficient of heat transfer";
    Modelica.Blocks.Interfaces.RealInput WindSpeedPort annotation(
      Placement(visible = true, transformation(origin = {-98, -72}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(extent = {{-2, -8}, {18, 12}}, rotation = 0)));
  protected
    Modelica.SIunits.Temp_C posDiff = abs(port_b.T - port_a.T) "Positive temperature difference";
  equation
// configure wind speed port
    if cardinality(WindSpeedPort) < 1 then
      WindSpeedPort = v;
    end if;
// no storage of heat
    port_a.Q_flow + port_b.Q_flow = 0;
/*
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      port_b -> wall
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      port_a -> air
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  */
// ------------------------------------------------------
// ceiling
    if control_type == hor_up then
      alpha = 0.5 * posDiff ^ 0.31;
// ------------------------------------------------------
// ground
    elseif control_type == hor_down then
      alpha = 2 * posDiff ^ 0.31;
//-------------------------------------------------------
// vertical interior wall
    elseif control_type == vert_interior then
      alpha = 1.6 * posDiff ^ 0.3;
//-------------------------------------------------------
// forced convection, outside the building
    elseif control_type == vert_exterior then
      alpha = 4 + 4 * WindSpeedPort;
// custom
//-------------------------------------------------------
    else
      alpha = alpha_custom;
    end if;
    port_a.Q_flow = alpha * A * (port_a.T - port_b.T);
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{0, 60}, {20, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {156, 156, 156}, fillPattern = FillPattern.Solid), Rectangle(extent = {{20, 60}, {40, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {182, 182, 182}, fillPattern = FillPattern.Solid), Rectangle(extent = {{40, 60}, {60, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {207, 207, 207}, fillPattern = FillPattern.Solid), Rectangle(extent = {{60, 60}, {80, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {244, 244, 244}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-80, 60}, {0, -100}}, lineColor = {0, 255, 255}, fillColor = {211, 243, 255}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Polygon(points = {{80, 60}, {80, 60}, {60, 20}, {60, 60}, {80, 60}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, lineThickness = 0.5, smooth = Smooth.None, fillColor = {157, 166, 208}, fillPattern = FillPattern.Solid), Polygon(points = {{60, 60}, {60, 20}, {40, -20}, {40, 60}, {60, 60}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, lineThickness = 0.5, smooth = Smooth.None, fillColor = {102, 110, 139}, fillPattern = FillPattern.Solid), Polygon(points = {{40, 60}, {40, -20}, {20, -60}, {20, 60}, {40, 60}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, lineThickness = 0.5, smooth = Smooth.None, fillColor = {75, 82, 103}, fillPattern = FillPattern.Solid), Polygon(points = {{20, 60}, {20, -60}, {0, -100}, {0, 60}, {20, 60}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, lineThickness = 0.5, smooth = Smooth.None, fillColor = {51, 56, 70}, fillPattern = FillPattern.Solid), Line(points = {{-20, 16}, {-20, -64}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None), Line(points = {{-20, 16}, {-30, 4}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None), Line(points = {{-38, 16}, {-48, 4}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None), Line(points = {{-54, 16}, {-64, 4}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None), Line(points = {{-38, 16}, {-38, -64}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None), Line(points = {{-54, 16}, {-54, -64}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.None)}),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(lineColor = {0, 255, 255}, fillColor = {211, 243, 255}, fillPattern = FillPattern.Solid, extent = {{-80, 60}, {0, -100}}), Rectangle(extent = {{-80, 60}, {80, -100}}), Rectangle(lineColor = {0, 0, 255}, fillColor = {244, 244, 244}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{60, 60}, {80, -100}}), Rectangle(lineColor = {0, 0, 255}, fillColor = {207, 207, 207}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{40, 60}, {60, -100}}), Rectangle(lineColor = {0, 0, 255}, fillColor = {182, 182, 182}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{20, 60}, {40, -100}}), Rectangle(lineColor = {0, 0, 255}, fillColor = {156, 156, 156}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{0, 60}, {20, -100}}), Polygon(lineColor = {0, 0, 255}, fillColor = {157, 166, 208}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{80, 60}, {80, 60}, {60, 20}, {60, 60}, {80, 60}}), Polygon(lineColor = {0, 0, 255}, fillColor = {102, 110, 139}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{60, 60}, {60, 20}, {40, -20}, {40, 60}, {60, 60}}), Polygon(lineColor = {0, 0, 255}, fillColor = {75, 82, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{40, 60}, {40, -20}, {20, -60}, {20, 60}, {40, 60}}), Polygon(lineColor = {0, 0, 255}, fillColor = {51, 56, 70}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{20, 60}, {20, -60}, {0, -100}, {0, 60}, {20, 60}}), Line(points = {{-58, 20}, {-68, 8}}, color = {0, 0, 255}, thickness = 0.5), Line(points = {{-58, 20}, {-58, -60}}, color = {0, 0, 255}, thickness = 0.5), Line(points = {{-40, 20}, {-50, 8}}, color = {0, 0, 255}, thickness = 0.5), Line(points = {{-40, 20}, {-40, -60}}, color = {0, 0, 255}, thickness = 0.5), Line(points = {{-22, 20}, {-32, 8}}, color = {0, 0, 255}, thickness = 0.5), Line(points = {{-22, 20}, {-22, -60}}, color = {0, 0, 255}, thickness = 0.5)}));
  end HeatConvUniversal;

  model GeneigteFlaeche
    import Modelica.SIunits.Conversions.from_deg;
    import Modelica.SIunits.Conversions.to_deg;
    parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg eps = 38.670 "latitude of location";
    parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg gamma = 13.4 "azimut of tilted surface, e.g. 0=south, 90=west, 180=north, -90=east";
    parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg beta = 90 "tilt of surface, e.g. 0=horizontal surface, 90=vertical surface";
    parameter Real rho_floor = 0.2 "ground reflection coefficient";
    Real cos_phi;
    Real phi;
    Real cos_phi_help;
    Real cos_theta_z;
    Real cos_theta_z_help;
    Real R;
    Real R_help;
    Real I_dir_hor "direct irradiation on horizontal surface";
    Real I_diff_hor "diffuse irradiation on horizontal surface";
    Real I_dir_phi "direct irradiation on tilted surface";
    Real I_diff_phi "diffuse irradiation on tilted surface";
    Real delta "declination";
    Real omega "hour angle";
    Modelica.Blocks.Interfaces.RealOutput Out_phi annotation(
      Placement(visible = true, transformation(origin = {90, -64}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.TwoRadiationPort out_rad annotation(
      Placement(visible = true, transformation(origin = {90, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.TwoRadiationPort in_rad annotation(
      Placement(visible = true, transformation(origin = {-90, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-96, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.SolarPositionPort solarposition annotation(
      Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    delta = solarposition.delta;
    omega = solarposition.omega;
// calculation of cos_theta_z [Duffie/Beckman, p.15], cos_theta_z is manually cut at 0 (no neg. values)
    cos_theta_z_help = sin(from_deg(delta)) * sin(from_deg(eps)) + cos(from_deg(delta)) * cos(from_deg(eps)) * cos(from_deg(omega));
    cos_theta_z = 0.5 * (cos_theta_z_help + abs(cos_theta_z_help));
// calculation of cos_phi [Duffie/Beckman, p.15], cos_phi is manually cut at 0 (no neg. values)
    cos_phi_help = sin(from_deg(delta)) * sin(from_deg(eps)) * cos(from_deg(beta)) - sin(from_deg(delta)) * cos(from_deg(eps)) * sin(from_deg(beta)) * cos(from_deg(gamma)) + cos(from_deg(delta)) * cos(from_deg(eps)) * cos(from_deg(beta)) * cos(from_deg(omega)) + cos(from_deg(delta)) * sin(from_deg(eps)) * sin(from_deg(beta)) * cos(from_deg(gamma)) * cos(from_deg(omega)) + cos(from_deg(delta)) * sin(from_deg(beta)) * sin(from_deg(gamma)) * sin(from_deg(omega));
    cos_phi = (cos_phi_help + abs(cos_phi_help)) / 2;
//cos_phi=cos_phi_help;
    phi = to_deg(acos(cos_phi));
    Out_phi = phi;
// calculation of R factor [Duffie/Beckman, p.25], due to numerical problems (cos_theta_z in denominator)
// R is manually set to 0 for theta_z >= 80° (-> 90° means sunset)
    if noEvent(cos_theta_z <= 0.17365) then
      R_help = cos_theta_z * cos_phi;
    else
      R_help = cos_phi / cos_theta_z;
    end if;
    R = R_help;
// calculation of total radiation on tilted surface according to model of Liu and Jordan
// according to [Dissertation Nytsch-Geusen, p.98]
    I_dir_hor = in_rad.I_dir;
    I_diff_hor = in_rad.I_diff;
    I_dir_phi = R * I_dir_hor;
    I_diff_phi = 0.5 * (1 + cos(from_deg(beta))) * I_diff_hor + 0.5 * rho_floor * (1 - cos(from_deg(beta))) * (I_dir_hor + I_diff_hor);
//Out_I_glob_phi.I = max(0, I_dir_phi + I_diff_phi);
    out_rad.I_diff = I_diff_phi;
    out_rad.I_dir = I_dir_phi;
    annotation(
      DymolaStoredErrors,
      Documentation(info = "<html>
<p><h4><font color=\"#008000\">Overview</font></h4></p>
<p>
The <b>RadOnTiltedSurf</b> model calculates the total radiance on a tilted surface.
</p>
<p><h4><font color=\"#008000\">Level of Development</font></h4></p>
<p><img src=\"modelica://HVAC/Images/stars4.png\"/></p>
<p><h4><font color=\"#008000\">Concept</font></h4></p>
<p>
The <b>RadOnTiltedSurf</b> model uses output data of the <a href=\"Sun\"><b>Sun</b></a> model and weather data (beam and diffuse radiance on a horizontal surface) to compute total radiance on a tilted surface. It needs information on the tilt angle and the azimut angle of the surface, the latitude of the location and the ground reflection coefficient.
</p>
<p><h4><font color=\"#008000\">Example Results</font></h4></p>
<p>The model is checked within the <a href=\"Building.Examples.Weather.WeatherModels\">weather</a> example as part of the <a href=\"Building.Components.Weather.Weather\">weather</a> model. </p>
</html>", revisions = "<html>
<ul>
  <li><i>May 02, 2013&nbsp;</i> by Ole Odendahl:<br/>Formatted documentation appropriately</li>
  <li><i>March 14, 2005&nbsp;</i>
         by Timo Haase:<br>
         Implemented.</li>
</ul>
</html>"),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {-6.79, -1}, extent = {{-44.16, 32.15}, {44.16, -32.15}}, textString = "geneigte Flaeche"), Rectangle(origin = {-1, 0}, extent = {{-85, 80}, {81, -100}})}),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-80, 60}, {80, -100}}), Rectangle(fillColor = {170, 213, 255}, pattern = LinePattern.None, fillPattern = FillPattern.HorizontalCylinder, extent = {{-80, 60}, {80, -100}}), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 225, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{14, 36}, {66, -16}}, endAngle = 360), Rectangle(fillColor = {0, 127, 0}, pattern = LinePattern.None, fillPattern = FillPattern.HorizontalCylinder, extent = {{-80, -40}, {80, -100}}), Rectangle(lineColor = {0, 0, 255}, fillColor = {0, 127, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-80, -72}, {80, -100}}), Polygon(fillColor = {226, 226, 226}, fillPattern = FillPattern.VerticalCylinder, points = {{-60, -64}, {-22, -76}, {-22, -32}, {-60, -24}, {-60, -64}}), Polygon(fillColor = {0, 77, 0}, pattern = LinePattern.None, fillPattern = FillPattern.VerticalCylinder, points = {{-60, -64}, {-80, -72}, {-80, -100}, {-60, -100}, {-22, -76}, {-60, -64}})}));
  end GeneigteFlaeche;

  model RadiativeHeatTransfer
    //import RefractionAngle;
    //import tau_a;
    //import tau;
    parameter Real N_pane = 1.526;
    parameter Integer number_panes = 2;
    parameter Real K = 4.0 "extinction coefficient";
    parameter Real d = 0.004;
    //parameter Real frameratio = 0.2;
    //parameter Real A = 2.0 "window area";
    Real I_dir_phi;
    Real I_diff_phi;
    Real tau_dir;
    Real tau_diff;
    Real tau_a_half_dir;
    Real tau_a_dir;
    Real tau_a_half_diff;
    Real tau_a_diff;
    Real phi_2;
    Real qdot_sw_in;
    Real qdot_sw_out;
    Real qdot_sw_a_in;
    Modelica.Blocks.Interfaces.RealInput In_phi_1 annotation(
      Placement(visible = true, transformation(origin = {-102, -52}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -52}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    BaseLib.TwoRadiationPort tworadiationport1 annotation(
      Placement(visible = true, transformation(origin = {-100, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.RadiationPort radiationport1 annotation(
      Placement(visible = true, transformation(origin = {98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    phi_2 = HomeAssignment_2023.RefractionAngle(In_phi_1);
    tau_dir = HomeAssignment_2023.tau(In_phi_1, d = d, K = K, number_panes = number_panes, N_pane = N_pane);
    tau_diff = HomeAssignment_2023.tau(60, d = d, K = K, number_panes = number_panes, N_pane = N_pane);
    tau_a_half_dir = HomeAssignment_2023.tau_a(In_phi_1, pos = 0.5, d = d, K = K, number_panes = number_panes, N_pane = N_pane);
    tau_a_dir = HomeAssignment_2023.tau_a(In_phi_1, d = d, K = K, number_panes = number_panes, N_pane = N_pane);
    tau_a_half_diff = HomeAssignment_2023.tau_a(60, pos = 0.5, d = d, K = K, number_panes = number_panes, N_pane = N_pane);
    tau_a_diff = HomeAssignment_2023.tau_a(60, d = d, K = K, number_panes = number_panes, N_pane = N_pane);
    tworadiationport1.I_dir = I_dir_phi;
    tworadiationport1.I_diff = I_diff_phi;
    qdot_sw_in = tau_dir * I_dir_phi + tau_diff * I_diff_phi;
    qdot_sw_out = qdot_sw_in * 0.1;
    qdot_sw_a_in = (tau_a_half_dir - tau_a_dir) * I_dir_phi + (tau_a_half_diff - tau_a_diff) * I_diff_phi + (1 - tau_a_half_diff) * qdot_sw_out;
    radiationport1.I = qdot_sw_in + qdot_sw_a_in;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end RadiativeHeatTransfer;

  function tau
    import Modelica.SIunits.Conversions.from_deg;
    input Real phi_1;
    input Integer number_panes = 2;
    input Real K = 4 "extinction coefficient";
    input Real d = 0.004 "pane thickness [m]";
    input Real N_pane;
    output Real tau_value;
  protected
    Real r_perpend;
    Real r_parallel;
    Real tau_r(max = 1);
    Real phi_2;
  algorithm
    phi_2 := HomeAssignment_2023.RefractionAngle(phi_1);
    r_perpend := sin(from_deg(phi_2 - phi_1)) ^ 2 / sin(from_deg(phi_2 + phi_1)) ^ 2;
    r_parallel := tan(from_deg(phi_2 - phi_1)) ^ 2 / tan(from_deg(phi_2 + phi_1)) ^ 2;
    tau_r := 0.5 * ((1. - r_parallel) / (1. + (2 * number_panes - 1) * r_parallel) + (1. - r_perpend) / (1. + (2 * number_panes - 1) * r_perpend));
    tau_value := HomeAssignment_2023.tau_a(phi_1, d = d, K = K, number_panes = number_panes, N_pane = N_pane) * tau_r;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end tau;

  function tau_a
    import Modelica.SIunits.Conversions.from_deg;
    import Modelica.SIunits.Conversions.to_deg;
    input Real phi_1;
    input Integer number_panes = 2;
    input Real K = 4 "extinction coefficient";
    input Real d = 0.004 "pane thickness [m]";
    input Real pos = 1 "spatial position in the pane [0,1]";
    input Real N_pane;
    output Real tau_a;
  protected
    Real phi_2;
  algorithm
    phi_2 := HomeAssignment_2023.RefractionAngle(phi_1, N_2 = N_pane);
    tau_a := exp(-pos * number_panes * d * K / cos(from_deg(phi_2)));
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end tau_a;

  function RefractionAngle
    import Modelica.SIunits.Conversions.to_deg;
    import Modelica.SIunits.Conversions.from_deg;
    input Real phi_1;
    input Real N_1 = 1.0;
    input Real N_2 = 1.526;
    output Real phi_2;
  algorithm
    if noEvent(from_deg(phi_1) <= 0.01) then
      phi_2 := to_deg(asin(0.01 * N_1 / N_2));
    else
      phi_2 := to_deg(asin(sin(from_deg(phi_1)) * N_1 / N_2));
    end if;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end RefractionAngle;

  model Wall
    parameter Real eps_LW "longwave emissivity inside";
    parameter Real alpha_SW = 0.7 "shortwave absorption outside";
    parameter Boolean outside = false "outside oder inside wall" annotation(
      Dialog(group = "relevant if vertical"));
    parameter Boolean withWindow = false "with window or not" annotation(
      Dialog(group = "relevant if vertical"));
    parameter Boolean roof = false "whether or not wall is a roof" annotation(
      Dialog(group = "relevant if horizontal"));
    parameter Boolean lowestFloor = false "whether or not wall is the lowest floor of the building" annotation(
      Dialog(group = "relevant if horizontal"));
    parameter Boolean horizontal "horizontal or vertical wall";
    parameter Modelica.SIunits.Height h = 3 "Height" annotation(
      Dialog(group = "Geometry"));
    parameter Modelica.SIunits.Length l = 4 "Length" annotation(
      Dialog(group = "Geometry"));
    parameter Integer n(min = 1) = 8 "Number of layers" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.Thickness d[n] = fill(0.007, n) "Thickness" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.Density rho[n] = fill(1600, n) "Density" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.ThermalConductivity lambda[n] = fill(2.4, n) "Thermal conductivity" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.SpecificHeatCapacity c[n] = fill(1000, n) "Specific heat capacity" annotation(
      Dialog(group = "Structure of wall layers", enable = not selectable));
    parameter Modelica.SIunits.Temperature T0 = Modelica.SIunits.Conversions.from_degC(16) "Initial temperature" annotation(
      Dialog(group = "Thermal"));
    parameter Real A_window = 2.0 "window area" annotation(
      Dialog(group = "Window"));
    parameter Real frameratio = 0.2 "window frame ratio" annotation(
      Dialog(group = "Window"));
    parameter Real N_pane = 1.526 "refraction index/Brechungsindex" annotation(
      Dialog(group = "Window"));
    parameter Real d_pane = 0.004 "pane thickness" annotation(
      Dialog(group = "Window"));
    parameter Integer number_panes = 2 "number of panes in the window" annotation(
      Dialog(group = "Window"));
    parameter Real beta = 90 "tilt angle" annotation(
      Dialog(group = "Orientation"));
    parameter Real gamma = 0 "azimuth angle" annotation(
      Dialog(group = "Orientation"));
    parameter Real rho_floor = 0.2 "ground reflection coefficient";
    parameter Real K = 4.0 "extinction parameter" annotation(
      Dialog(group = "Window"));
    parameter Real U = 2.1 "U value" annotation(
      Dialog(group = "Window"));
    parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg lat = 43.273 "latitude of location";
    Modelica.Blocks.Interfaces.RealInput u if (outside or roof) and not lowestFloor annotation(
      Placement(visible = true, transformation(origin = {-102, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-98, -52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.Star star2 if not outside and not lowestFloor and not roof annotation(
      Placement(visible = true, transformation(origin = {-100, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.window window1(A = A_window, frameratio = frameratio, N_pane = N_pane, d = d_pane, number_panes = number_panes, beta = beta, gamma = gamma, rho_floor = rho_floor, K = K, U = U, eps = lat) if outside and withWindow annotation(
      Placement(visible = true, transformation(origin = {1, 69}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
    BaseLib.Star star1 annotation(
      Placement(visible = true, transformation(origin = {98, 82}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.TwoRadiationPort tworadiationport1 if (outside or roof) and not lowestFloor annotation(
      Placement(visible = true, transformation(origin = {-100, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-98, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.SolarPositionPort solarposition if (outside or roof) and not lowestFloor annotation(
      Placement(visible = true, transformation(origin = {-102, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-96, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(visible = true, transformation(origin = {-98, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-98, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation(
      Placement(visible = true, transformation(origin = {98, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    parameter HomeAssignment_2023.HeatConvUniversal.Temp control_type_up_or_out = if not horizontal then if outside then HomeAssignment_2023.HeatConvUniversal.vert_exterior else HomeAssignment_2023.HeatConvUniversal.vert_interior else HomeAssignment_2023.HeatConvUniversal.hor_down "control type of outer/upper convection node";
    parameter HomeAssignment_2023.HeatConvUniversal.Temp control_type_down_or_in = if not horizontal then HomeAssignment_2023.HeatConvUniversal.vert_interior else HomeAssignment_2023.HeatConvUniversal.hor_up "control type of inner/lower convection node";
    BaseLib.addRadiation addradiation1(alpha_SW = alpha_SW) if (outside or roof) and not lowestFloor annotation(
      Placement(visible = true, transformation(origin = {-49, 59}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    HomeAssignment_2023.GeneigteFlaeche geneigteflaeche1(eps = lat, gamma = gamma, beta = beta, rho_floor = rho_floor) if (outside or roof) and not lowestFloor annotation(
      Placement(visible = true, transformation(origin = {-71, 55}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Blocks.Interfaces.IntegerInput WindowOpeningInput if withWindow "Integer input, 1 if window open, 0 if closed" annotation(
      Placement(transformation(extent = {{-112, 80}, {-92, 100}}), iconTransformation(extent = {{-112, 80}, {-92, 100}})));
  protected
    parameter Modelica.SIunits.Area A = h * l - A_window;
    HomeAssignment_2023.N_Layer n_layer1(T0 = T0, c = c, clearance = A_window, d = d, h = h, l = l, lambda = lambda, n = n, rho = rho) annotation(
      Placement(visible = true, transformation(origin = {2, 4}, extent = {{-36, -36}, {36, 36}}, rotation = 0)));
    HomeAssignment_2023.HeatConvUniversal heatconvuniversal1(control_type = control_type_up_or_out, A = A) annotation(
      Placement(visible = true, transformation(origin = {-64, 4}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    BaseLib.TwoStar_RadEx twostar_radex1(A = A, eps = eps_LW) annotation(
      Placement(visible = true, transformation(origin = {48, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.TwoStar_RadEx twostar_radex2(A = A, eps = 0.9) if not outside and not lowestFloor and not roof annotation(
      Placement(visible = true, transformation(origin = {-58, -66}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.HeatConvUniversal heatconvuniversal2(control_type = control_type_up_or_out, A = A) annotation(
      Placement(visible = true, transformation(origin = {66, 4}, extent = {{14, -14}, {-14, 14}}, rotation = 0)));
    BaseLib.I_to_Qdot i_to_qdot1(A = A) if (outside or roof) and not lowestFloor annotation(
      Placement(visible = true, transformation(origin = {-34, 60}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  equation
    connect(geneigteflaeche1.out_rad, addradiation1.tworadiationport1) annotation(
      Line(points = {{-61.76, 57.42}, {-58, 57.42}, {-58, 58}, {-55.86, 58}, {-55.86, 59}}));
    connect(solarposition, geneigteflaeche1.solarposition) annotation(
      Line(points = {{-102, 70}, {-92, 70}, {-92, 58.52}, {-80.9, 58.52}}, color = {0, 0, 127}));
    connect(tworadiationport1, geneigteflaeche1.in_rad) annotation(
      Line(points = {{-100, 42}, {-86, 42}, {-86, 47.3}, {-81.56, 47.3}}, color = {255, 128, 0}));
    connect(addradiation1.radiationport1, i_to_qdot1.radiationport1) annotation(
      Line(points = {{-42, 59}, {-39.58, 59}, {-39.58, 60.06}}));
    connect(i_to_qdot1.port_a, n_layer1.port_a) annotation(
      Line(points = {{-28.6, 59.88}, {-34, 59.88}, {-34, 4.72}}, color = {191, 0, 0}));
    connect(window1.star1, star1) annotation(
      Line(points = {{16.3, 78.3}, {52, 78.3}, {52, 82}, {98, 82}, {98, 82}}, color = {95, 95, 95}));
    connect(tworadiationport1, window1.tworadiationport1) annotation(
      Line(points = {{-100, 42}, {-86, 42}, {-86, 34}, {-26, 34}, {-26, 72}, {-14.9, 72}, {-14.9, 71.7}}, color = {255, 128, 0}));
    connect(solarposition, window1.solarposition) annotation(
      Line(points = {{-102, 70}, {-40, 70}, {-40, 78}, {-13.1, 78}, {-13.1, 78.9}}, color = {0, 0, 127}));
    connect(window1.heatport_inside, port_b) annotation(
      Line(points = {{16, 57.9}, {26, 57.9}, {26, 20}, {86, 20}, {86, 4}, {98, 4}}, color = {191, 0, 0}));
    connect(heatconvuniversal2.port_a, port_b) annotation(
      Line(points = {{80, 4.28}, {90, 4.28}, {90, 4}, {98, 4}}, color = {191, 0, 0}));
    connect(n_layer1.port_b, heatconvuniversal2.port_b) annotation(
      Line(points = {{38, 4}, {52, 4}}, color = {191, 0, 0}));
    connect(port_a, window1.heatport_outside) annotation(
      Line(points = {{-98, 4}, {-84, 4}, {-84, 20}, {-22, 20}, {-22, 58}, {-14.3, 58}, {-14.3, 57.9}}, color = {191, 0, 0}));
    connect(heatconvuniversal1.port_a, port_a) annotation(
      Line(points = {{-76, 4.24}, {-88, 4.24}, {-88, 4}, {-98, 4}}, color = {191, 0, 0}));
    connect(twostar_radex1.star, star1) annotation(
      Line(points = {{57.8, 60}, {80, 60}, {80, 82}, {98, 82}}, color = {95, 95, 95}));
    connect(n_layer1.port_b, twostar_radex1.port_a) annotation(
      Line(points = {{38, 4}, {38, 60.4}}, color = {191, 0, 0}));
    connect(twostar_radex2.port_a, n_layer1.port_a) annotation(
      Line(points = {{-48, -65.6}, {-34, -65.6}, {-34, 4.72}, {-34, 4.72}}, color = {191, 0, 0}));
    connect(star2, twostar_radex2.star) annotation(
      Line(points = {{-100, -90}, {-68, -90}, {-68, -66}, {-67.8, -66}}, color = {0, 0, 255}));
    connect(u, heatconvuniversal1.WindSpeedPort) annotation(
      Line(points = {{-102, -30}, {-80, -30}, {-80, 4.24}, {-63.04, 4.24}}, color = {0, 0, 127}));
    connect(heatconvuniversal1.port_b, n_layer1.port_a) annotation(
      Line(points = {{-52, 4}, {-44, 4}, {-44, 4.72}, {-34, 4.72}}, color = {191, 0, 0}));
    connect(window1.u, WindowOpeningInput) annotation(
      Line(points = {{-13.7, 66}, {-52, 66}, {-52, 90}, {-102, 90}}, color = {255, 127, 0}));
    connect(star2, star2) annotation(
      Line(points = {{-100, -90}, {-100, -90}, {-100, -90}}, color = {95, 95, 95}));
    connect(port_a, port_a) annotation(
      Line(points = {{-98, 4}, {-98, 4}}, color = {191, 0, 0}));
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {1, 0}, fillColor = {149, 149, 149}, fillPattern = FillPattern.Solid, extent = {{-45, 92}, {41, -92}}), Rectangle(origin = {-45, 0}, fillColor = {202, 228, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-51, 92}, {1, -92}}), Rectangle(origin = {93, 0}, fillColor = {202, 228, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-51, 92}, {1, -92}}), Line(origin = {-72, 64}, points = {{-14, 12}, {18, -8}, {18, -8}}, color = {255, 255, 0}, thickness = 4), Line(origin = {-64, 80}, points = {{-18, 8}, {18, -8}, {18, -8}}, color = {255, 255, 0}, thickness = 4), Line(origin = {-78, 44}, points = {{-10, 20}, {18, -8}, {18, -8}}, color = {255, 255, 0}, thickness = 4)}));
  end Wall;

  model SolarPositionModel
    // This models calculates the solar declination delta and hour angle omega as a function of longitude and time.
    parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg L = 6.640 "longitude of location in degree";
    Real day "day of the year";
    Real t_sol "solar time";
    Real Z "correction term for solar time";
    Real x "auxiliary variable";
    Real delta "solar declination";
    Real omega "hour angle";
    Real omega_help "hour angle (auxiliary)";
    constant Real pi = Modelica.Constants.pi;
    BaseLib.SolarPositionPort solarposition annotation(
      Placement(visible = true, transformation(origin = {-10, -18}, extent = {{80, 10}, {100, 30}}, rotation = 0), iconTransformation(origin = {103, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  equation
    day = integer(time / 86400) + 1;
    x = 360 / 365.26 * day - 2.72;
    if time < 7430400 or time > 26179200 then
      (t_sol - time) / 3600.0 = (-4 / 60.0 * (15 - L)) + 1 / 60.0 * Z;
    else
      (t_sol - time) / 3600.0 = (-4 / 60.0 * (30 - L)) + 1 / 60.0 * Z;
    end if;
    Z = (-7.66 * sin(x * pi / 180)) - 9.87 * sin((2 * x + 24.99 + 3.83 * sin(x * pi / 180)) * pi / 180);
    omega_help = 15 * (t_sol - 12 * 3600) / 3600;
    omega = mod(omega_help, 360);
    delta = 23.45 * sin((284 + day) / 365.26 * (pi / 180) * 360);
    solarposition.delta = delta;
    solarposition.omega = omega;
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0, grid = {2, 2})),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0, grid = {2, 2}), graphics = {Ellipse(origin = {-3, -10}, lineColor = {255, 255, 0}, fillColor = {255, 180, 0}, fillPattern = FillPattern.Solid, extent = {{-53, 54}, {53, -54}}, endAngle = 360), Ellipse(origin = {13, 7}, fillPattern = FillPattern.Solid, extent = {{-11, 23}, {11, -23}}, endAngle = 360), Line(origin = {-11, -36}, points = {{-25, 6}, {-11, -4}}, thickness = 2), Ellipse(origin = {-17, 7}, fillPattern = FillPattern.Solid, extent = {{-11, 23}, {11, -23}}, endAngle = 360), Line(origin = {-75, 37}, points = {{-21, 3}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {-57, 55}, points = {{-15, 15}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {-79, 7}, points = {{-15, -13}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {-77, -21}, points = {{-17, -25}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {-33, -59}, points = {{1, -37}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {-27, 67}, points = {{3, 29}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {11, 63}, points = {{31, 31}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {47, 21}, points = {{49, -11}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {45, 49}, points = {{51, 23}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {47, -35}, points = {{45, -33}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {-61, -45}, points = {{-19, -41}, {15, -15}, {15, -15}}, color = {255, 180, 0}, thickness = 3), Line(origin = {45, -75}, points = {{-7, 7}, {15, -15}, {21, -21}}, color = {255, 180, 0}, thickness = 3), Line(origin = {3, -46}, points = {{-25, 6}, {-9, 2}}, thickness = 2), Line(origin = {19, -50}, points = {{-25, 6}, {-11, 8}}, thickness = 2), Line(origin = {33, -48}, points = {{-25, 6}, {-13, 10}}, thickness = 2), Line(origin = {45, -44}, points = {{-25, 6}, {-17, 12}}, thickness = 2), Line(origin = {-15, -38}, points = {{-25, 6}, {-17, 12}}, thickness = 2), Line(origin = {55, -40}, points = {{-25, 6}, {-31, 14}}, thickness = 2)}));
  end SolarPositionModel;

  class OrthoRad "Orthogonal radiation exchange"
    extends BaseLib.TwoPort;
    parameter Modelica.SIunits.Length edge = 4 "common edge";
    parameter Modelica.SIunits.Length side1 = 4 "side of wall 1";
    parameter Modelica.SIunits.Length side2 = 4 "side of wall 2";
    parameter Modelica.SIunits.Emissivity eps1 = 0.95 "emission coefficient wall 1";
    parameter Modelica.SIunits.Emissivity eps2 = 0.95 "emission coefficient wall 2";
    constant Real pi = Modelica.Constants.pi "declaration of pi";
    Real F12;
  protected
    parameter Real B = side1 / edge;
    parameter Real C = side2 / edge;
    parameter Real B2 = B ^ 2;
    parameter Real C2 = C ^ 2;
    parameter Real A1 = side1 * edge;
    parameter Real A2 = side2 * edge;
  equation
// view factors according to VDI Wärmeatlas:
    F12 = 1 / (pi * B) * (B * Modelica.Math.atan(1 / B) + C * Modelica.Math.atan(1 / C) - sqrt(B2 + C2) * Modelica.Math.atan(1 / sqrt(B2 + C2)) + 1 / 4 * (B2 * Modelica.Math.log((1 + B2 + C2) * B2 / ((B2 + C2) * (1 + B2))) + C2 * Modelica.Math.log((1 + B2 + C2) * C2 / ((B2 + C2) * (1 + C2))) - Modelica.Math.log((1 + B2 + C2) / ((1 + B2) * (1 + C2)))));
    port_a.Q_flow = F12 * eps1 * eps2 * A1 * Modelica.Constants.sigma * (Modelica.SIunits.Conversions.from_degC(port_a.T) ^ 4 - Modelica.SIunits.Conversions.from_degC(port_b.T) ^ 4);
    port_b.Q_flow = -port_a.Q_flow;
    annotation(
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Rectangle(extent = {{-80, 60}, {60, 40}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {191, 95, 0}), Rectangle(extent = {{60, 40}, {80, -100}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.VerticalCylinder, fillColor = {191, 95, 0}), Rectangle(extent = {{-80, 40}, {60, -100}}, lineColor = {0, 0, 0}, fillColor = {159, 223, 223}, fillPattern = FillPattern.Solid), Rectangle(extent = {{60, 60}, {80, 40}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-48, -6}, {4, -62}}, color = {255, 255, 0}), Line(points = {{-34, 4}, {18, -52}}, color = {255, 255, 0}), Polygon(points = {{18, -52}, {12, -38}, {6, -46}, {18, -52}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid), Polygon(points = {{-50, -4}, {-46, -16}, {-40, -10}, {-50, -4}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid)}),
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {60, 40}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {191, 95, 0}), Rectangle(extent = {{60, 40}, {80, -100}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.VerticalCylinder, fillColor = {191, 95, 0}), Rectangle(extent = {{-80, 40}, {60, -100}}, lineColor = {0, 0, 0}, fillColor = {159, 223, 223}, fillPattern = FillPattern.Solid), Rectangle(extent = {{60, 60}, {80, 40}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-48, -6}, {4, -62}}, color = {255, 255, 0}), Line(points = {{-34, 4}, {18, -52}}, color = {255, 255, 0}), Polygon(points = {{18, -52}, {12, -38}, {6, -46}, {18, -52}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid), Polygon(points = {{-50, -4}, {-46, -16}, {-40, -10}, {-50, -4}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid)}),
      Window(x = 0.07000000000000001, y = 0.03, width = 0.99, height = 0.65));
  end OrthoRad;

  class ParallelRad "Parallel radiation exchange"
    extends BaseLib.TwoPort;
    parameter Modelica.SIunits.Length side1 = 4 "common width of walls";
    parameter Modelica.SIunits.Length side2 = 4 "common length of walls";
    parameter Modelica.SIunits.Length distance = 4 "distance between walls";
    parameter Modelica.SIunits.Emissivity eps1 = 0.95 "emission coefficient wall 1";
    parameter Modelica.SIunits.Emissivity eps2 = 0.95 "emission coefficient wall 2";
    constant Real pi = Modelica.Constants.pi "declaration of pi";
    Real Phi12;
  protected
    parameter Real B = side1 / distance;
    parameter Real C = side2 / distance;
    parameter Real B2 = B ^ 2;
    parameter Real BC = B * C;
    parameter Real C2 = C ^ 2;
    parameter Real A1 = side1 * side2;
    parameter Real A2 = side1 * side2;
  equation
// View factors according to VDI Wärmeatlas:
    Phi12 = 1 / pi * (1 / BC * log((1 + B2) * (1 + C2) / (1 + B2 + C2)) - 2 / B * Modelica.Math.atan(C) - 2 / C * Modelica.Math.atan(B) + 2 / C * sqrt(1 + C2) * Modelica.Math.atan(B / sqrt(1 + C2)) + 2 / B * sqrt(1 + B2) * Modelica.Math.atan(C / sqrt(1 + B2)));
    port_a.Q_flow = Phi12 * eps1 * eps2 * A1 * Modelica.Constants.sigma * (Modelica.SIunits.Conversions.from_degC(port_a.T) ^ 4 - Modelica.SIunits.Conversions.from_degC(port_b.T) ^ 4);
    port_b.Q_flow = -port_a.Q_flow;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {-60, -100}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.VerticalCylinder, fillColor = {191, 95, 0}), Rectangle(extent = {{-60, -100}, {60, 60}}, lineColor = {0, 0, 0}, fillColor = {159, 223, 223}, fillPattern = FillPattern.Solid), Line(points = {{-40, 0}, {40, 0}}, color = {255, 255, 0}), Rectangle(extent = {{60, 60}, {80, -100}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.VerticalCylinder, fillColor = {191, 95, 0}), Line(points = {{-40, -40}, {40, -40}}, color = {255, 255, 0}), Polygon(points = {{-40, 0}, {-30, 6}, {-30, -6}, {-40, 0}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid), Line(points = {{-40, 20}, {40, 20}}, color = {255, 255, 0}), Polygon(points = {{40, -40}, {32, -34}, {32, -46}, {40, -40}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}), Rectangle(extent = {{-80, 60}, {-60, -100}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.VerticalCylinder, fillColor = {191, 95, 0}), Rectangle(extent = {{-60, -100}, {60, 60}}, lineColor = {0, 0, 0}, fillColor = {159, 223, 223}, fillPattern = FillPattern.Solid), Line(points = {{-40, 0}, {40, 0}}, color = {255, 255, 0}), Rectangle(extent = {{60, 60}, {80, -100}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.VerticalCylinder, fillColor = {191, 95, 0}), Line(points = {{-40, -40}, {40, -40}}, color = {255, 255, 0}), Polygon(points = {{-38, 0}, {-28, 6}, {-28, -6}, {-38, 0}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid), Line(points = {{-40, 20}, {40, 20}}, color = {255, 255, 0}), Polygon(points = {{40, -40}, {32, -34}, {32, -46}, {40, -40}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid)}),
      DymolaStoredErrors);
  end ParallelRad;

  model RadExchange
    parameter Modelica.SIunits.Length length = 4 "length of walls 1 and 3";
    parameter Modelica.SIunits.Length width = 4 "width of walls 2 and 4";
    parameter Modelica.SIunits.Length height = 4 "height of room";
    parameter Modelica.SIunits.Emissivity epsNorth = 0.95 "emission coefficient wall 1";
    parameter Modelica.SIunits.Emissivity epsEast = 0.95 "emission coefficient wall 2";
    parameter Modelica.SIunits.Emissivity epsSouth = 0.95 "emission coefficient wall 3";
    parameter Modelica.SIunits.Emissivity epsWest = 0.95 "emission coefficient wall 4";
    parameter Modelica.SIunits.Emissivity epsC = 0.95 "emission coefficient ceiling";
    parameter Modelica.SIunits.Emissivity epsF = 0.95 "emission coefficient floor";
    HomeAssignment_2023.OrthoRad north_west(edge = height, side1 = width, side2 = length, eps1 = epsNorth, eps2 = epsWest) annotation(
      Placement(visible = true, transformation(origin = {-60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad north_east(edge = height, side1 = width, side2 = length, eps1 = epsNorth, eps2 = epsEast) annotation(
      Placement(visible = true, transformation(origin = {-20, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad ceiling_north(edge = width, side1 = length, side2 = height, eps1 = epsC, eps2 = epsNorth) annotation(
      Placement(visible = true, transformation(origin = {40, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad ceiling_west(edge = length, side1 = width, side2 = height, eps1 = epsC, eps2 = epsWest) annotation(
      Placement(visible = true, transformation(origin = {20, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad floor_north(edge = width, side1 = length, side2 = height, eps1 = epsF, eps2 = epsNorth) annotation(
      Placement(visible = true, transformation(origin = {40, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad floor_west(edge = length, side1 = width, side2 = height, eps1 = epsF, eps2 = epsWest) annotation(
      Placement(visible = true, transformation(origin = {20, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad south_east(edge = height, side1 = width, side2 = length, eps1 = epsSouth, eps2 = epsEast) annotation(
      Placement(visible = true, transformation(origin = {-20, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad south_west(edge = height, side1 = width, side2 = length, eps1 = epsSouth, eps2 = epsWest) annotation(
      Placement(visible = true, transformation(origin = {-60, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad floor_east(edge = length, side1 = width, side2 = height, eps1 = epsF, eps2 = epsEast) annotation(
      Placement(visible = true, transformation(origin = {60, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad ceiling_east(edge = length, side1 = width, side2 = height, eps1 = epsC, eps2 = epsEast) annotation(
      Placement(visible = true, transformation(origin = {60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad floor_south(edge = width, side1 = length, side2 = height, eps1 = epsF, eps2 = epsSouth) annotation(
      Placement(visible = true, transformation(origin = {80, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.OrthoRad ceiling_south(edge = width, side1 = length, side2 = height, eps1 = epsC, eps2 = epsSouth) annotation(
      Placement(visible = true, transformation(origin = {80, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.ParallelRad north_south(side1 = width, side2 = height, distance = length, eps1 = epsNorth, eps2 = epsSouth) annotation(
      Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a west annotation(
      Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a north annotation(
      Placement(visible = true, transformation(origin = {-40, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a south annotation(
      Placement(visible = true, transformation(origin = {-40, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a ceiling annotation(
      Placement(visible = true, transformation(origin = {60, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a floor annotation(
      Placement(visible = true, transformation(origin = {60, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.ParallelRad ceiling_floor(side1 = width, side2 = length, distance = height, eps1 = epsC, eps2 = epsF) annotation(
      Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    HomeAssignment_2023.ParallelRad west_east(side1 = length, side2 = height, distance = width, eps1 = epsWest, eps2 = epsEast) annotation(
      Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a east annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ceiling_north.port_b, north) annotation(
      Line(color = {255, 128, 0}, points = {{40, 30}, {40, 25.7951}, {4.24028, 25.7951}, {4.24028, 78.7986}, {-40, 78.7986}, {-40, 80}}));
    connect(ceiling_south.port_b, south) annotation(
      Line(color = {255, 128, 0}, points = {{80, 30}, {80, 30.0353}, {95.7597, 30.0353}, {95.7597, -95.053}, {-39.9293, -95.053}, {-39.9293, -80}, {-40, -80}}));
    connect(floor_south.port_a, south) annotation(
      Line(color = {255, 128, 0}, points = {{80.2, -30}, {80.2, -24.3816}, {95.7597, -24.3816}, {95.7597, -95.053}, {-39.9293, -95.053}, {-39.9293, -80}, {-40, -80}}));
    connect(floor_east.port_a, east) annotation(
      Line(color = {255, 128, 0}, points = {{60.2, -30}, {60.2, -18.0212}, {29.3286, -18.0212}, {29.3286, 1.06007}, {0, 1.06007}, {0, 0}}));
    connect(ceiling_east.port_b, east) annotation(
      Line(color = {255, 128, 0}, points = {{60, 30}, {60, 15.9011}, {29.682, 15.9011}, {29.682, 1.06007}, {0, 1.06007}, {0, 0}}));
    connect(floor_north.port_a, north) annotation(
      Line(color = {255, 128, 0}, points = {{40.2, -30}, {40.2, 21.9081}, {4.24028, 21.9081}, {4.24028, 78.7986}, {-43.8163, 78.7986}, {-43.8163, 80}, {-40, 80}}));
    connect(floor_west.port_a, west) annotation(
      Line(color = {255, 128, 0}, points = {{20.2, -30}, {20.2, 25.7951}, {-79.5053, 25.7951}, {-79.5053, 4.947}, {-80, 4.947}, {-80, 0}}));
    connect(ceiling_west.port_b, west) annotation(
      Line(color = {255, 128, 0}, points = {{20, 30}, {20, 25.7951}, {-79.5053, 25.7951}, {-79.5053, 0}, {-80, 0}}));
    connect(north_east.port_b, east) annotation(
      Line(points = {{-20, 30}, {-20, 24.0283}, {-0.706714, 24.0283}, {-0.706714, 0}, {0, 0}}));
    connect(south, south_east.port_b) annotation(
      Line(points = {{-40, -80}, {-39.9293, -80}, {-39.9293, -54.7703}, {-19.4346, -54.7703}, {-19.4346, -50}, {-20, -50}}));
    connect(south, south_west.port_b) annotation(
      Line(points = {{-40, -80}, {-39.9293, -80}, {-39.9293, -54.417}, {-59.364, -54.417}, {-59.364, -50}, {-60, -50}}));
    connect(south, north_south.port_b) annotation(
      Line(points = {{-40, -80}, {-39.9293, -80}, {-39.9293, -10}, {-40, -10}}));
    connect(east, south_east.port_a) annotation(
      Line(points = {{0, 0}, {-1.41343, 0}, {-1.41343, -14.1343}, {-19.4346, -14.1343}, {-19.4346, -30}, {-19.8, -30}}));
    connect(east, west_east.port_b) annotation(
      Line(points = {{0, 0}, {-1.41343, 0}, {-1.41343, -14.1343}, {-19.4346, -14.1343}, {-19.4346, -10}, {-20, -10}}));
    connect(west, south_west.port_a) annotation(
      Line(points = {{-80, 0}, {-61.1307, 0}, {-61.1307, -30.3887}, {-59.8, -30.3887}, {-59.8, -30}}));
    connect(west, west_east.port_a) annotation(
      Line(points = {{-80, 0}, {-61.1307, 0}, {-61.1307, 15.9011}, {-20.1413, 15.9011}, {-20.1413, 9.89399}, {-19.8, 9.89399}, {-19.8, 10}}));
    connect(west, north_west.port_b) annotation(
      Line(points = {{-80, 0}, {-61.1307, 0}, {-61.1307, 29.3286}, {-60, 29.3286}, {-60, 30}}));
    connect(north, north_east.port_a) annotation(
      Line(points = {{-40, 80}, {-39.576, 80}, {-39.576, 54.0636}, {-19.788, 54.0636}, {-19.788, 50}, {-19.8, 50}}));
    connect(north, north_west.port_a) annotation(
      Line(points = {{-40, 80}, {-39.576, 80}, {-39.576, 54.417}, {-60.424, 54.417}, {-60.424, 50}, {-59.8, 50}}));
    connect(north, north_south.port_a) annotation(
      Line(points = {{-40, 80}, {-39.576, 80}, {-39.576, 10}, {-39.8, 10}}));
    connect(ceiling_floor.port_a, ceiling) annotation(
      Line(points = {{60.2, 10}, {60.2, 10.2473}, {92.5795, 10.2473}, {92.5795, 53.3569}, {60.0707, 53.3569}, {60.0707, 80}, {60, 80}}));
    connect(ceiling_floor.port_b, floor) annotation(
      Line(points = {{60, -10}, {60, -10.2473}, {92.2261, -10.2473}, {92.2261, -54.7703}, {60.0707, -54.7703}, {60.0707, -80}, {60, -80}}));
    connect(floor, floor_west.port_b) annotation(
      Line(points = {{60, -80}, {60.0707, -80}, {60.0707, -54.7703}, {20.1413, -54.7703}, {20.1413, -50}, {20, -50}}));
    connect(floor, floor_north.port_b) annotation(
      Line(points = {{60, -80}, {60.0707, -80}, {60.0707, -54.7703}, {39.9293, -54.7703}, {39.9293, -50}, {40, -50}}));
    connect(floor, floor_south.port_b) annotation(
      Line(points = {{60, -80}, {60.0707, -80}, {60.0707, -54.7703}, {80.212, -54.7703}, {80.212, -50}, {80, -50}}));
    connect(floor, floor_east.port_b) annotation(
      Line(points = {{60, -80}, {60.0707, -80}, {60.0707, -50}, {60, -50}}));
    connect(ceiling, ceiling_west.port_a) annotation(
      Line(points = {{60, 80}, {60.0707, 80}, {60.0707, 53.3569}, {20.1413, 53.3569}, {20.1413, 50}, {20.2, 50}}));
    connect(ceiling, ceiling_north.port_a) annotation(
      Line(points = {{60, 80}, {60.0707, 80}, {60.0707, 53.3569}, {40.2827, 53.3569}, {40.2827, 50}, {40.2, 50}}));
    connect(ceiling, ceiling_south.port_a) annotation(
      Line(points = {{60, 80}, {60.0707, 80}, {60.0707, 53.3569}, {80.212, 53.3569}, {80.212, 50}, {80.2, 50}}));
    connect(ceiling, ceiling_east.port_a) annotation(
      Line(points = {{60, 80}, {60.0707, 80}, {60.0707, 50}, {60.2, 50}}));
    annotation(
      Documentation(revisions = "<html>
<ul>
  <li><i>March 14, 2005&nbsp;</i>
         by Timo Haase:<br>
         Implemented.</li>
</ul>
</html>", info = "<html>
<p>
The <b>RadExchange</b> model describes energy transfer by radiation in a typical room with four walls (simple box model). There are three pairs of identical surfaces: wall 1/wall 3, wall 2/wall 4 and ceiling/floor.
</p>

<dl>
<dt><b>Main Author:</b>
<dd>Timo Haase <br>
    Technische Universtit&auml;t Berlin <br>
    Hermann-Rietschel-Institut <br>
    Marchstr. 4 <br> 
    D-10587 Berlin <br>
    e-mail: <a href=\"mailto:timo.haase@tu-berlin.de\">timo.haase@tu-berlin.de</a><br>
</dl>

</html>"),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Polygon(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, points = {{-80, 40}, {-40, 80}, {40, 80}, {80, 40}, {80, -40}, {40, -80}, {-40, -80}, {-80, -40}, {-80, 40}}), Line(points = {{-52, 56}, {38, 56}, {-52, 56}, {-52, -36}, {-52, 56}, {0, -52}, {-52, 56}, {62, -26}, {-52, 56}}), Line(points = {{-44, -56}, {48, -56}, {-44, -56}, {40, 30}, {-44, -56}, {-4, 46}, {-44, -56}, {-44, 18}}, color = {255, 0, 0}), Text(origin = {-82.86, 22.61}, extent = {{-14.66, 8.48}, {23.85, -16.61}}, textString = "west"), Text(origin = {-75.23, 93.85}, extent = {{-14.66, 8.48}, {23.85, -16.61}}, textString = "north"), Text(origin = {-73.46, -76.47}, extent = {{-14.66, 8.48}, {23.85, -16.61}}, textString = "south"), Text(origin = {70.57, 26.57}, extent = {{-14.66, 8.48}, {23.85, -16.61}}, textString = "east"), Text(origin = {54.53, 84.03}, extent = {{-14.66, 8.48}, {23.85, -16.61}}, textString = "ceiling"), Text(origin = {54.3936, -76.8866}, extent = {{-14.66, 8.48}, {23.85, -16.61}}, textString = "floor")}));
  end RadExchange;

  class VarHeatFlow "Variable heat flow source"
    extends BaseLib.OnePort;
    parameter Modelica.SIunits.HeatFlowRate Q_flow = 2000 "heat flow rate";
    Modelica.Blocks.Interfaces.RealInput InPort annotation(
      Placement(transformation(extent = {{-100, -90}, {-80, -70}}, rotation = 0)));
  equation
    port_a.Q_flow = -Q_flow * InPort;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}, fillColor = {159, 191, 223}, fillPattern = FillPattern.Solid), Polygon(points = {{-12, -46}, {-22, -32}, {-26, -18}, {-24, -6}, {-8, 10}, {0, 16}, {10, 32}, {20, 14}, {24, -8}, {22, -24}, {16, -38}, {10, -46}, {-12, -46}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid), Polygon(points = {{-8, -46}, {-14, -34}, {-16, -22}, {-12, -12}, {0, -2}, {10, 12}, {14, 0}, {18, -14}, {14, -30}, {8, -46}, {-8, -46}}, lineColor = {255, 127, 0}, fillColor = {255, 127, 0}, fillPattern = FillPattern.Solid), Polygon(points = {{-4, -46}, {-10, -32}, {-8, -20}, {2, -12}, {6, -6}, {10, -20}, {10, -30}, {8, -36}, {4, -46}, {-4, -46}}, lineColor = {255, 0, 0}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid), Polygon(points = {{-4, -46}, {-6, -34}, {-4, -24}, {2, -18}, {4, -30}, {2, -40}, {0, -46}, {-4, -46}}, lineColor = {127, 0, 255}, fillColor = {0, 0, 191}, fillPattern = FillPattern.Solid), Polygon(points = {{-40, -46}, {40, -46}, {24, -70}, {-24, -70}, {-40, -46}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-80, 60}, {80, -100}}, lineColor = {0, 0, 0}, fillColor = {159, 191, 223}, fillPattern = FillPattern.Solid), Polygon(points = {{-12, -46}, {-22, -32}, {-26, -18}, {-24, -6}, {-8, 10}, {0, 16}, {10, 32}, {20, 14}, {24, -8}, {22, -24}, {16, -38}, {10, -46}, {-12, -46}}, lineColor = {255, 255, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid), Polygon(points = {{-8, -46}, {-14, -34}, {-16, -22}, {-12, -12}, {0, -2}, {10, 12}, {14, 0}, {18, -14}, {14, -30}, {8, -46}, {-8, -46}}, lineColor = {255, 127, 0}, fillColor = {255, 127, 0}, fillPattern = FillPattern.Solid), Polygon(points = {{-4, -46}, {-10, -32}, {-8, -20}, {2, -12}, {6, -6}, {10, -20}, {10, -30}, {8, -36}, {4, -46}, {-4, -46}}, lineColor = {255, 0, 0}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid), Polygon(points = {{-4, -46}, {-6, -34}, {-4, -24}, {2, -18}, {4, -30}, {2, -40}, {0, -46}, {-4, -46}}, lineColor = {127, 0, 255}, fillColor = {0, 0, 191}, fillPattern = FillPattern.Solid), Polygon(points = {{-40, -46}, {40, -46}, {24, -70}, {-24, -70}, {-40, -46}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid)}),
      Window(x = 0.23, y = 0.12, width = 0.6, height = 0.6),
      Documentation(revisions = "<html>
Taken from ATplus library.
</html>", info = "<html>
<p>
The <b>VarHeatFlow</b> model represents a variable heat flow source, controlled by an input port. The given heat flow rate is being multiplied by the signal on the input port and the weight factor.
</p>

</html>"));
  end VarHeatFlow;

  model ImportWeatherModified
    parameter String weatherfile = "C:\CR_Protyping Project DSS - BIM\Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype\openModelica\Saint-Tropez-min.dat";
    String line "Variable for each line of the text file";
    Boolean endOfFile(start = false) "True, if end of file is reached";
    Integer nextIndex(start = 1) "Index for the horizontal position in a line";
    Integer month "month";
    Integer dm "day in month (1-31)";
    Integer dy "day in year (1-365)";
    Integer h "hour (1-24)";
    Integer min "minute (1-60)";
    Real I_glob "global irradiance [W]";
    Real S_height "sun height [°]";
    Real G_gex "Extraterrestrial radiation [W]";
    Real I_glob_hr "global horizontal irradiance (considers a high horizon) [W]";
    Real I_diff "diffuse irradiance [W]";
    Real I_glob_inclined "inclined global irradiance (2pi) [W]";
    Real I_diff_inclined "inclined diffused irradiance (2pi) [W]";
    Real I_dir_beam "direct irradiance [W]";
    Real T_a "ambient temperature";
    Real FF "wind speed [m/s]";
    Real I_glob_bc "global irradiance with clear sky[W]";
    Real T_b "dew point temperature";
    final constant Real pi = 2 * Modelica.Math.asin(1.0);
    // 3.14159265358979;
    final constant Real D2R = pi / 180 "Degree to Radian";
    BaseLib.TwoRadiationPort tworadiationport1 annotation(
      Placement(visible = true, transformation(origin = {98, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(visible = true, transformation(origin = {102, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput windspeed annotation(
      Placement(visible = true, transformation(origin = {102, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.IntegerOutput HourOfDay "from 0 to 23" annotation(
      Placement(transformation(extent = {{94, 56}, {114, 76}})));
  algorithm
// For time steps 0+i*60, i=0,1,... do:
    when sample(0, 60) then
      (line, endOfFile) := Modelica.Utilities.Streams.readLine(weatherfile, integer(time / 60) + 4);
      (month, nextIndex) := Modelica.Utilities.Strings.scanInteger(line);
      (dm, nextIndex) := Modelica.Utilities.Strings.scanInteger(line, nextIndex + 1);
      (dy, nextIndex) := Modelica.Utilities.Strings.scanInteger(line, nextIndex + 1);
      (h, nextIndex) := Modelica.Utilities.Strings.scanInteger(line, nextIndex + 1);
      (min, nextIndex) := Modelica.Utilities.Strings.scanInteger(line, nextIndex + 1);
      (I_glob, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (S_height, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (G_gex, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (I_glob_hr, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (I_diff, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (I_glob_inclined, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (I_diff_inclined, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (I_dir_beam, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (T_a, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (FF, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (I_glob_bc, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
      (T_b, nextIndex) := Modelica.Utilities.Strings.scanReal(line, nextIndex + 1);
    end when;
  equation
    tworadiationport1.I_diff = I_diff;
    tworadiationport1.I_dir = I_dir_beam * sin(D2R * S_height);
    HourOfDay = h;
    port_a.T = T_a + 273.15;
    windspeed = FF;
    connect(HourOfDay, HourOfDay) annotation(
      Line(points = {{104, 66}, {104, 66}}, color = {255, 127, 0}));
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {9, -14}, extent = {{-89, 66}, {65, -52}}, textString = "read weather data"), Rectangle(origin = {2, -14}, extent = {{-96, 82}, {96, -82}})}),
      experiment(StartTime = 0, StopTime = 3.1536e+07, Tolerance = 1e-06, Interval = 10));
  end ImportWeatherModified;

  model Window_Status_hour_of_day
    Modelica.Blocks.Interfaces.IntegerOutput Window_Status annotation(
      Placement(visible = true, transformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.IntegerInput HourOfDay annotation(
      Placement(transformation(extent = {{76, 26}, {98, 48}})));
  equation
    if HourOfDay >= 6 and HourOfDay < 17 then
      Window_Status = 1;
// Closed
    elseif HourOfDay >= 18 and HourOfDay < 21 then
      Window_Status = 2;
// Tilted
    elseif HourOfDay >= 22 or HourOfDay < 5 then
      Window_Status = 3;
// Fully Open
    else
      Window_Status = 1;
// Closed (default)
    end if;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Window_Status_hour_of_day;

  model PrototypeSampleCase_SaintTropez
    HomeAssignment_2023.Wall wall_North(A_window = 2.25, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = 180, h = 3.2, horizontal = false, l = 4, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_East(A_window = 9, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = -90, h = 3.2, horizontal = false, l = 6, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_South(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 0, h = 3.2, horizontal = false, l = 4, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_west(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 90, h = 3.2, horizontal = false, l = 6, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall roof(alpha_SW = 0, beta = 0, c = {950, 900, 1100}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.25, 0.10, 0.01}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {2.0, 0.043, 0.415}, lat = 43.273, lowestFloor = false, n = 3, outside = true, rho(each displayUnit = "kg/m3") = {2400, 100, 1248}, rho_floor = 0.2, roof = true, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {30, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall floor(T0(displayUnit = "K") = 289.15, alpha_SW = 0, beta = 0, c = {1600, 900, 950}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.021, 0.12, 0.30}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {0.13, 0.043, 2.0}, lat = 43.273, lowestFloor = true, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {462, 100, 2400}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.SolarPositionModel solarPositionModel(L = 43.273) annotation(
      Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow = 0.01, T_ref = 289.15, alpha = 1) annotation(
      Placement(visible = true, transformation(origin = {-90, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.Airload airload(T(fixed = true), T0(displayUnit = "K") = 289.15, V = 76.8, c = 1007, rho = 1.19) annotation(
      Placement(visible = true, transformation(origin = {-10, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(visible = true, transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-12, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ImportWeatherModified importWeatherModified annotation(
      Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature fixedTemperature1(T = 16) annotation(
      Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.RadExchange radExchange(epsC = 0.95, epsEast = 0.95, epsF = 0.95, epsNorth = 0.95, epsSouth = 0.95, epsWest = 0.95, height = 3.2, length = 6, width = 4) annotation(
      Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Window_Status_hour_of_day window_Status_hour_of_day annotation(
      Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(fixedHeatFlow.port, floor.port_a) annotation(
      Line(points = {{-80, -50}, {-60, -50}, {-60, -50}, {-60, -50}}, color = {191, 0, 0}));
    connect(solarPositionModel.solarposition, roof.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, 36}, {40, 36}, {40, 36}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_North.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -4}, {40, -4}, {40, -4}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_East.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -44}, {40, -44}, {40, -44}}, color = {0, 0, 127}));
    connect(wall_South.port_b, port_a) annotation(
      Line(points = {{-40, 30}, {-34, 30}, {-34, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, floor.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -50}, {-40, -50}}, color = {191, 0, 0}));
    connect(airload.port_a, port_a) annotation(
      Line(points = {{-10, -60}, {-10, -60}, {-10, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, wall_west.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -10}, {-40, -10}}, color = {191, 0, 0}));
    connect(port_a, wall_East.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -50}, {20, -50}, {20, -50}, {20, -50}}, color = {191, 0, 0}));
    connect(port_a, wall_North.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -10}, {20, -10}, {20, -10}}, color = {191, 0, 0}));
    connect(port_a, roof.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, 30}, {20, 30}, {20, 30}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, wall_East.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -56}, {40, -56}, {40, -56}}, color = {0, 0, 127}));
    connect(importWeatherModified.windspeed, wall_North.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -16}, {40, -16}, {40, -16}}, color = {0, 0, 127}));
    connect(importWeatherModified.port_a, wall_North.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -12}, {40, -12}, {40, -10}}, color = {191, 0, 0}));
    connect(importWeatherModified.tworadiationport1, wall_East.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -48}, {40, -48}, {40, -48}}));
    connect(importWeatherModified.tworadiationport1, wall_North.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -8}, {40, -8}, {40, -8}}));
    connect(importWeatherModified.port_a, wall_East.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -52}, {40, -52}, {40, -50}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, roof.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, 24}, {40, 24}, {40, 24}}, color = {0, 0, 127}));
    connect(importWeatherModified.tworadiationport1, roof.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, 32}, {40, 32}, {40, 32}}));
    connect(importWeatherModified.port_a, roof.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, 28}, {40, 28}, {40, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_South.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, 30}, {-60, 30}, {-60, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_west.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -10}, {-60, -10}, {-60, -10}}, color = {191, 0, 0}));
    connect(radExchange.north, wall_North.star1) annotation(
      Line(points = {{-14, 20}, {-14, 36}, {6, 36}, {6, -4}, {20, -4}}, color = {191, 0, 0}));
    connect(radExchange.west, wall_west.star1) annotation(
      Line(points = {{-20, 10}, {-28, 10}, {-28, -4}, {-40, -4}, {-40, -4}}, color = {191, 0, 0}));
    connect(radExchange.floor, floor.star1) annotation(
      Line(points = {{-4, 0}, {-4, 0}, {-4, -10}, {-28, -10}, {-28, -44}, {-40, -44}, {-40, -44}}, color = {191, 0, 0}));
    connect(radExchange.ceiling, roof.star1) annotation(
      Line(points = {{-4, 20}, {6, 20}, {6, 36}, {20, 36}, {20, 36}, {20, 36}}, color = {191, 0, 0}));
    connect(radExchange.east, wall_East.star1) annotation(
      Line(points = {{0, 10}, {6, 10}, {6, -44}, {20, -44}, {20, -44}}, color = {191, 0, 0}));
    connect(radExchange.south, wall_South.star1) annotation(
      Line(points = {{-14, 0}, {-14, 0}, {-14, -4}, {-28, -4}, {-28, 36}, {-40, 36}, {-40, 36}}, color = {191, 0, 0}));
    connect(window_Status_hour_of_day.HourOfDay, importWeatherModified.HourOfDay) annotation(
      Line(points = {{38, 74}, {56, 74}, {56, 56}, {76, 56}, {76, 60}, {76, 60}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day.Window_Status, wall_North.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -2}, {40, -2}, {40, 0}, {40, 0}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day.Window_Status, wall_East.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -42}, {40, -42}, {40, -40}}, color = {255, 127, 0}));
  protected
    annotation(
      experiment(StartTime = 15638400, StopTime = 16848000, Tolerance = 0.001, Interval = 10));
  end PrototypeSampleCase_SaintTropez;

  model Window_Status_hour_of_day_Winter
    Modelica.Blocks.Interfaces.IntegerOutput Window_Status annotation(
      Placement(visible = true, transformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.IntegerInput HourOfDay annotation(
      Placement(transformation(extent = {{76, 26}, {98, 48}})));
  equation
    if HourOfDay >= 10 and HourOfDay < 13 then
      Window_Status = 2;
// Tilted
//elseif HourOfDay >= 18 and HourOfDay < 21 then
//Window_Status = 2; // Tilted
//elseif HourOfDay >= 22 or HourOfDay < 5 then
//Window_Status = 3; // Fully Open
    else
      Window_Status = 1;
// Closed (default)
    end if;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Window_Status_hour_of_day_Winter;

  model Window_Status_hour_of_day_Spring
    Modelica.Blocks.Interfaces.IntegerOutput Window_Status annotation(
      Placement(visible = true, transformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.IntegerInput HourOfDay annotation(
      Placement(transformation(extent = {{76, 26}, {98, 48}})));
  equation
    if HourOfDay >= 10 and HourOfDay < 13 then
      Window_Status = 2;
// Closed
    elseif HourOfDay >= 14 and HourOfDay < 16 then
      Window_Status = 3;
// Tilted
//elseif HourOfDay >= 22 or HourOfDay < 5 then
//Window_Status = 3; // Fully Open
    else
      Window_Status = 1;
// Closed (default)
    end if;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Window_Status_hour_of_day_Spring;

  model Window_Status_hour_of_day_Summer
    Modelica.Blocks.Interfaces.IntegerOutput Window_Status annotation(
      Placement(visible = true, transformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.IntegerInput HourOfDay annotation(
      Placement(transformation(extent = {{76, 26}, {98, 48}})));
  equation
    if HourOfDay >= 8 and HourOfDay < 12 then
      Window_Status = 3;
// Fully Open
    elseif HourOfDay >= 13 and HourOfDay < 20 then
      Window_Status = 2;
// Tilted
//elseif HourOfDay >= 22 or HourOfDay < 5 then
//Window_Status = 3; // Fully Open
    else
      Window_Status = 1;
// Closed (default)
    end if;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Window_Status_hour_of_day_Summer;

  model Window_Status_hour_of_day_Autumn
    Modelica.Blocks.Interfaces.IntegerOutput Window_Status annotation(
      Placement(visible = true, transformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.IntegerInput HourOfDay annotation(
      Placement(transformation(extent = {{76, 26}, {98, 48}})));
  equation
    if HourOfDay >= 8 and HourOfDay < 10 then
      Window_Status = 3;
// Closed
    elseif HourOfDay >= 11 and HourOfDay < 17 then
      Window_Status = 2;
// Tilted
//elseif HourOfDay >= 22 or HourOfDay < 5 then
//Window_Status = 3; // Fully Open
    else
      Window_Status = 1;
// Closed (default)
    end if;
    annotation(
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Window_Status_hour_of_day_Autumn;

  model PrototypeSampleCase_SaintTropez_Spring
    HomeAssignment_2023.Wall wall_North(A_window = 2.25, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = 180, h = 3.2, horizontal = false, l = 4, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_East(A_window = 9, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = -90, h = 3.2, horizontal = false, l = 6, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_South(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 0, h = 3.2, horizontal = false, l = 4, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_west(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 90, h = 3.2, horizontal = false, l = 6, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall roof(alpha_SW = 0.95, beta = 0, c = {950, 900, 1100}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.25, 0.10, 0.01}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {2.0, 0.043, 0.415}, lat = 43.273, lowestFloor = false, n = 3, outside = true, rho(each displayUnit = "kg/m3") = {2400, 100, 1248}, rho_floor = 0.2, roof = true, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {30, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall floor(T0(displayUnit = "K") = 289.15, alpha_SW = 0, beta = 0, c = {1600, 900, 950}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.021, 0.12, 0.30}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {0.13, 0.043, 2.0}, lat = 43.273, lowestFloor = true, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {462, 100, 2400}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.SolarPositionModel solarPositionModel(L = 43.273) annotation(
      Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    BaseLib.Airload airload(T(fixed = true), T0(displayUnit = "K") = 289.15, V = 76.8, c = 1007, rho = 1.19) annotation(
      Placement(visible = true, transformation(origin = {-10, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(visible = true, transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-12, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ImportWeatherModified importWeatherModified annotation(
      Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature fixedTemperature1(T = 16) annotation(
      Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.RadExchange radExchange(epsC = 0.95, epsEast = 0.95, epsF = 0.95, epsNorth = 0.95, epsSouth = 0.95, epsWest = 0.95, height = 3.2, length = 6, width = 4) annotation(
      Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Window_Status_hour_of_day_Spring window_Status_hour_of_day_Spring annotation(
      Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(solarPositionModel.solarposition, roof.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, 36}, {40, 36}, {40, 36}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_North.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -4}, {40, -4}, {40, -4}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_East.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -44}, {40, -44}, {40, -44}}, color = {0, 0, 127}));
    connect(wall_South.port_b, port_a) annotation(
      Line(points = {{-40, 30}, {-34, 30}, {-34, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, floor.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -50}, {-40, -50}}, color = {191, 0, 0}));
    connect(airload.port_a, port_a) annotation(
      Line(points = {{-10, -60}, {-10, -60}, {-10, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, wall_west.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -10}, {-40, -10}}, color = {191, 0, 0}));
    connect(port_a, wall_East.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -50}, {20, -50}, {20, -50}, {20, -50}}, color = {191, 0, 0}));
    connect(port_a, wall_North.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -10}, {20, -10}, {20, -10}}, color = {191, 0, 0}));
    connect(port_a, roof.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, 30}, {20, 30}, {20, 30}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, wall_East.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -56}, {40, -56}, {40, -56}}, color = {0, 0, 127}));
    connect(importWeatherModified.windspeed, wall_North.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -16}, {40, -16}, {40, -16}}, color = {0, 0, 127}));
    connect(importWeatherModified.port_a, wall_North.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -12}, {40, -12}, {40, -10}}, color = {191, 0, 0}));
    connect(importWeatherModified.tworadiationport1, wall_East.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -48}, {40, -48}, {40, -48}}));
    connect(importWeatherModified.tworadiationport1, wall_North.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -8}, {40, -8}, {40, -8}}));
    connect(importWeatherModified.port_a, wall_East.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -52}, {40, -52}, {40, -50}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, roof.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, 24}, {40, 24}, {40, 24}}, color = {0, 0, 127}));
    connect(importWeatherModified.tworadiationport1, roof.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, 32}, {40, 32}, {40, 32}}));
    connect(importWeatherModified.port_a, roof.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, 28}, {40, 28}, {40, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_South.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, 30}, {-60, 30}, {-60, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_west.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -10}, {-60, -10}, {-60, -10}}, color = {191, 0, 0}));
    connect(radExchange.north, wall_North.star1) annotation(
      Line(points = {{-14, 20}, {-14, 36}, {6, 36}, {6, -4}, {20, -4}}, color = {191, 0, 0}));
    connect(radExchange.west, wall_west.star1) annotation(
      Line(points = {{-20, 10}, {-28, 10}, {-28, -4}, {-40, -4}, {-40, -4}}, color = {191, 0, 0}));
    connect(radExchange.floor, floor.star1) annotation(
      Line(points = {{-4, 0}, {-4, 0}, {-4, -10}, {-28, -10}, {-28, -44}, {-40, -44}, {-40, -44}}, color = {191, 0, 0}));
    connect(radExchange.ceiling, roof.star1) annotation(
      Line(points = {{-4, 20}, {6, 20}, {6, 36}, {20, 36}, {20, 36}, {20, 36}}, color = {191, 0, 0}));
    connect(radExchange.east, wall_East.star1) annotation(
      Line(points = {{0, 10}, {6, 10}, {6, -44}, {20, -44}, {20, -44}}, color = {191, 0, 0}));
    connect(radExchange.south, wall_South.star1) annotation(
      Line(points = {{-14, 0}, {-14, 0}, {-14, -4}, {-28, -4}, {-28, 36}, {-40, 36}, {-40, 36}}, color = {191, 0, 0}));
    connect(window_Status_hour_of_day_Spring.HourOfDay, importWeatherModified.HourOfDay) annotation(
      Line(points = {{38, 74}, {56, 74}, {56, 56}, {76, 56}, {76, 60}, {76, 60}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Spring.Window_Status, wall_North.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -2}, {40, -2}, {40, 0}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Spring.Window_Status, wall_East.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -42}, {40, -42}, {40, -40}}, color = {255, 127, 0}));
    connect(fixedTemperature1.port, floor.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -50}, {-60, -50}, {-60, -50}}, color = {191, 0, 0}));
  protected
    annotation(
      experiment(StartTime = 5.0976e+06, StopTime = 1.30464e+07, Tolerance = 0.001, Interval = 10));
  end PrototypeSampleCase_SaintTropez_Spring;

  model PrototypeSampleCase_SaintTropez_Summer
    HomeAssignment_2023.Wall wall_North(A_window = 2.25, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = 180, h = 3.2, horizontal = false, l = 4, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_East(A_window = 9, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = -90, h = 3.2, horizontal = false, l = 6, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_South(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 0, h = 3.2, horizontal = false, l = 4, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_west(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 90, h = 3.2, horizontal = false, l = 6, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall roof(alpha_SW = 0.95, beta = 0, c = {950, 900, 1100}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.25, 0.10, 0.01}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {2.0, 0.043, 0.415}, lat = 43.273, lowestFloor = false, n = 3, outside = true, rho(each displayUnit = "kg/m3") = {2400, 100, 1248}, rho_floor = 0.2, roof = true, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {30, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall floor(T0(displayUnit = "K") = 289.15, alpha_SW = 0, beta = 0, c = {1600, 900, 950}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.021, 0.12, 0.30}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {0.13, 0.043, 2.0}, lat = 43.273, lowestFloor = true, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {462, 100, 2400}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.SolarPositionModel solarPositionModel(L = 43.273) annotation(
      Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    BaseLib.Airload airload(T(fixed = true), T0(displayUnit = "K") = 289.15, V = 76.8, c = 1007, rho = 1.19) annotation(
      Placement(visible = true, transformation(origin = {-10, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(visible = true, transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-12, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ImportWeatherModified importWeatherModified annotation(
      Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature fixedTemperature1(T = 16) annotation(
      Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.RadExchange radExchange(epsC = 0.95, epsEast = 0.95, epsF = 0.95, epsNorth = 0.95, epsSouth = 0.95, epsWest = 0.95, height = 3.2, length = 6, width = 4) annotation(
      Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Window_Status_hour_of_day_Summer window_Status_hour_of_day_Summer annotation(
      Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(solarPositionModel.solarposition, roof.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, 36}, {40, 36}, {40, 36}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_North.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -4}, {40, -4}, {40, -4}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_East.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -44}, {40, -44}, {40, -44}}, color = {0, 0, 127}));
    connect(wall_South.port_b, port_a) annotation(
      Line(points = {{-40, 30}, {-34, 30}, {-34, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, floor.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -50}, {-40, -50}}, color = {191, 0, 0}));
    connect(airload.port_a, port_a) annotation(
      Line(points = {{-10, -60}, {-10, -60}, {-10, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, wall_west.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -10}, {-40, -10}}, color = {191, 0, 0}));
    connect(port_a, wall_East.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -50}, {20, -50}, {20, -50}, {20, -50}}, color = {191, 0, 0}));
    connect(port_a, wall_North.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -10}, {20, -10}, {20, -10}}, color = {191, 0, 0}));
    connect(port_a, roof.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, 30}, {20, 30}, {20, 30}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, wall_East.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -56}, {40, -56}, {40, -56}}, color = {0, 0, 127}));
    connect(importWeatherModified.windspeed, wall_North.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -16}, {40, -16}, {40, -16}}, color = {0, 0, 127}));
    connect(importWeatherModified.port_a, wall_North.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -12}, {40, -12}, {40, -10}}, color = {191, 0, 0}));
    connect(importWeatherModified.tworadiationport1, wall_East.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -48}, {40, -48}, {40, -48}}));
    connect(importWeatherModified.tworadiationport1, wall_North.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -8}, {40, -8}, {40, -8}}));
    connect(importWeatherModified.port_a, wall_East.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -52}, {40, -52}, {40, -50}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, roof.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, 24}, {40, 24}, {40, 24}}, color = {0, 0, 127}));
    connect(importWeatherModified.tworadiationport1, roof.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, 32}, {40, 32}, {40, 32}}));
    connect(importWeatherModified.port_a, roof.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, 28}, {40, 28}, {40, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_South.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, 30}, {-60, 30}, {-60, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_west.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -10}, {-60, -10}, {-60, -10}}, color = {191, 0, 0}));
    connect(radExchange.north, wall_North.star1) annotation(
      Line(points = {{-14, 20}, {-14, 36}, {6, 36}, {6, -4}, {20, -4}}, color = {191, 0, 0}));
    connect(radExchange.west, wall_west.star1) annotation(
      Line(points = {{-20, 10}, {-28, 10}, {-28, -4}, {-40, -4}, {-40, -4}}, color = {191, 0, 0}));
    connect(radExchange.floor, floor.star1) annotation(
      Line(points = {{-4, 0}, {-4, 0}, {-4, -10}, {-28, -10}, {-28, -44}, {-40, -44}, {-40, -44}}, color = {191, 0, 0}));
    connect(radExchange.ceiling, roof.star1) annotation(
      Line(points = {{-4, 20}, {6, 20}, {6, 36}, {20, 36}, {20, 36}, {20, 36}}, color = {191, 0, 0}));
    connect(radExchange.east, wall_East.star1) annotation(
      Line(points = {{0, 10}, {6, 10}, {6, -44}, {20, -44}, {20, -44}}, color = {191, 0, 0}));
    connect(radExchange.south, wall_South.star1) annotation(
      Line(points = {{-14, 0}, {-14, 0}, {-14, -4}, {-28, -4}, {-28, 36}, {-40, 36}, {-40, 36}}, color = {191, 0, 0}));
    connect(window_Status_hour_of_day_Summer.HourOfDay, importWeatherModified.HourOfDay) annotation(
      Line(points = {{38, 74}, {56, 74}, {56, 56}, {76, 56}, {76, 60}, {76, 60}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Summer.Window_Status, wall_North.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, 0}, {40, 0}, {40, 0}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Summer.Window_Status, wall_East.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -42}, {40, -42}, {40, -40}}, color = {255, 127, 0}));
    connect(fixedTemperature1.port, floor.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -50}, {-60, -50}, {-60, -50}}, color = {191, 0, 0}));
  protected
    annotation(
      experiment(StartTime = 1.30464e+07, StopTime = 2.09952e+07, Tolerance = 0.001, Interval = 10));
  end PrototypeSampleCase_SaintTropez_Summer;

  model PrototypeSampleCase_SaintTropez_Autumn
    HomeAssignment_2023.Wall wall_North(A_window = 2.25, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = 180, h = 3.2, horizontal = false, l = 4, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_East(A_window = 9, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = -90, h = 3.2, horizontal = false, l = 6, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_South(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 0, h = 3.2, horizontal = false, l = 4, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_west(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 90, h = 3.2, horizontal = false, l = 6, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall roof(alpha_SW = 0.95, beta = 0, c = {950, 900, 1100}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.25, 0.10, 0.01}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {2.0, 0.043, 0.415}, lat = 43.273, lowestFloor = false, n = 3, outside = true, rho(each displayUnit = "kg/m3") = {2400, 100, 1248}, rho_floor = 0.2, roof = true, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {30, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall floor(T0(displayUnit = "K") = 289.15, alpha_SW = 0, beta = 0, c = {1600, 900, 950}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.021, 0.12, 0.30}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {0.13, 0.043, 2.0}, lat = 43.273, lowestFloor = true, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {462, 100, 2400}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.SolarPositionModel solarPositionModel(L = 43.273) annotation(
      Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    BaseLib.Airload airload(T(fixed = true), T0(displayUnit = "K") = 289.15, V = 76.8, c = 1007, rho = 1.19) annotation(
      Placement(visible = true, transformation(origin = {-10, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(visible = true, transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-12, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ImportWeatherModified importWeatherModified annotation(
      Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature fixedTemperature1(T = 16) annotation(
      Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.RadExchange radExchange(epsC = 0.95, epsEast = 0.95, epsF = 0.95, epsNorth = 0.95, epsSouth = 0.95, epsWest = 0.95, height = 3.2, length = 6, width = 4) annotation(
      Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Window_Status_hour_of_day_Autumn window_Status_hour_of_day_Autumn annotation(
      Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(solarPositionModel.solarposition, roof.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, 36}, {40, 36}, {40, 36}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_North.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -4}, {40, -4}, {40, -4}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_East.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -44}, {40, -44}, {40, -44}}, color = {0, 0, 127}));
    connect(wall_South.port_b, port_a) annotation(
      Line(points = {{-40, 30}, {-34, 30}, {-34, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, floor.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -50}, {-40, -50}}, color = {191, 0, 0}));
    connect(airload.port_a, port_a) annotation(
      Line(points = {{-10, -60}, {-10, -60}, {-10, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, wall_west.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -10}, {-40, -10}}, color = {191, 0, 0}));
    connect(port_a, wall_East.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -50}, {20, -50}, {20, -50}, {20, -50}}, color = {191, 0, 0}));
    connect(port_a, wall_North.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -10}, {20, -10}, {20, -10}}, color = {191, 0, 0}));
    connect(port_a, roof.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, 30}, {20, 30}, {20, 30}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, wall_East.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -56}, {40, -56}, {40, -56}}, color = {0, 0, 127}));
    connect(importWeatherModified.windspeed, wall_North.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -16}, {40, -16}, {40, -16}}, color = {0, 0, 127}));
    connect(importWeatherModified.port_a, wall_North.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -12}, {40, -12}, {40, -10}}, color = {191, 0, 0}));
    connect(importWeatherModified.tworadiationport1, wall_East.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -48}, {40, -48}, {40, -48}}));
    connect(importWeatherModified.tworadiationport1, wall_North.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -8}, {40, -8}, {40, -8}}));
    connect(importWeatherModified.port_a, wall_East.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -52}, {40, -52}, {40, -50}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, roof.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, 24}, {40, 24}, {40, 24}}, color = {0, 0, 127}));
    connect(importWeatherModified.tworadiationport1, roof.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, 32}, {40, 32}, {40, 32}}));
    connect(importWeatherModified.port_a, roof.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, 28}, {40, 28}, {40, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_South.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, 30}, {-60, 30}, {-60, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_west.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -10}, {-60, -10}, {-60, -10}}, color = {191, 0, 0}));
    connect(radExchange.north, wall_North.star1) annotation(
      Line(points = {{-14, 20}, {-14, 36}, {6, 36}, {6, -4}, {20, -4}}, color = {191, 0, 0}));
    connect(radExchange.west, wall_west.star1) annotation(
      Line(points = {{-20, 10}, {-28, 10}, {-28, -4}, {-40, -4}, {-40, -4}}, color = {191, 0, 0}));
    connect(radExchange.floor, floor.star1) annotation(
      Line(points = {{-4, 0}, {-4, 0}, {-4, -10}, {-28, -10}, {-28, -44}, {-40, -44}, {-40, -44}}, color = {191, 0, 0}));
    connect(radExchange.ceiling, roof.star1) annotation(
      Line(points = {{-4, 20}, {6, 20}, {6, 36}, {20, 36}, {20, 36}, {20, 36}}, color = {191, 0, 0}));
    connect(radExchange.east, wall_East.star1) annotation(
      Line(points = {{0, 10}, {6, 10}, {6, -44}, {20, -44}, {20, -44}}, color = {191, 0, 0}));
    connect(radExchange.south, wall_South.star1) annotation(
      Line(points = {{-14, 0}, {-14, 0}, {-14, -4}, {-28, -4}, {-28, 36}, {-40, 36}, {-40, 36}}, color = {191, 0, 0}));
    connect(window_Status_hour_of_day_Autumn.HourOfDay, importWeatherModified.HourOfDay) annotation(
      Line(points = {{38, 74}, {56, 74}, {56, 54}, {76, 54}, {76, 60}, {76, 60}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Autumn.Window_Status, wall_North.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -2}, {40, -2}, {40, 0}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Autumn.Window_Status, wall_East.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -42}, {40, -42}, {40, -40}}, color = {255, 127, 0}));
    connect(fixedTemperature1.port, floor.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -52}, {-60, -52}, {-60, -50}}, color = {191, 0, 0}));
  protected
    annotation(
      experiment(StartTime = 2.09952e+07, StopTime = 2.88576e+07, Tolerance = 0.001, Interval = 10));
  end PrototypeSampleCase_SaintTropez_Autumn;

  model PrototypeSampleCase_SaintTropez_WinterJanFeb
    HomeAssignment_2023.Wall wall_North(A_window = 2.25, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = 180, h = 3.2, horizontal = false, l = 4, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_East(A_window = 9, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = -90, h = 3.2, horizontal = false, l = 6, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_South(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 0, h = 3.2, horizontal = false, l = 4, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_west(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 90, h = 3.2, horizontal = false, l = 6, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall roof(alpha_SW = 0.95, beta = 0, c = {950, 900, 1100}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.25, 0.10, 0.01}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {2.0, 0.043, 0.415}, lat = 43.273, lowestFloor = false, n = 3, outside = true, rho(each displayUnit = "kg/m3") = {2400, 100, 1248}, rho_floor = 0.2, roof = true, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {30, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall floor(T0(displayUnit = "K") = 289.15, alpha_SW = 0, beta = 0, c = {1600, 900, 950}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.021, 0.12, 0.30}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {0.13, 0.043, 2.0}, lat = 43.273, lowestFloor = true, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {462, 100, 2400}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.SolarPositionModel solarPositionModel(L = 43.273) annotation(
      Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow = 0.01, T_ref = 289.15, alpha = 1) annotation(
      Placement(visible = true, transformation(origin = {-90, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.Airload airload(T(fixed = true), T0(displayUnit = "K") = 289.15, V = 76.8, c = 1007, rho = 1.19) annotation(
      Placement(visible = true, transformation(origin = {-10, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(visible = true, transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-12, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ImportWeatherModified importWeatherModified annotation(
      Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature fixedTemperature1(T = 16) annotation(
      Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.RadExchange radExchange(epsC = 0.95, epsEast = 0.95, epsF = 0.95, epsNorth = 0.95, epsSouth = 0.95, epsWest = 0.95, height = 3.2, length = 6, width = 4) annotation(
      Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Window_Status_hour_of_day_Winter window_Status_hour_of_day_Winter annotation(
      Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(fixedHeatFlow.port, floor.port_a) annotation(
      Line(points = {{-80, -50}, {-60, -50}, {-60, -50}, {-60, -50}}, color = {191, 0, 0}));
    connect(solarPositionModel.solarposition, roof.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, 36}, {40, 36}, {40, 36}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_North.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -4}, {40, -4}, {40, -4}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_East.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -44}, {40, -44}, {40, -44}}, color = {0, 0, 127}));
    connect(wall_South.port_b, port_a) annotation(
      Line(points = {{-40, 30}, {-34, 30}, {-34, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, floor.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -50}, {-40, -50}}, color = {191, 0, 0}));
    connect(airload.port_a, port_a) annotation(
      Line(points = {{-10, -60}, {-10, -60}, {-10, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, wall_west.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -10}, {-40, -10}}, color = {191, 0, 0}));
    connect(port_a, wall_East.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -50}, {20, -50}, {20, -50}, {20, -50}}, color = {191, 0, 0}));
    connect(port_a, wall_North.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -10}, {20, -10}, {20, -10}}, color = {191, 0, 0}));
    connect(port_a, roof.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, 30}, {20, 30}, {20, 30}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, wall_East.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -56}, {40, -56}, {40, -56}}, color = {0, 0, 127}));
    connect(importWeatherModified.windspeed, wall_North.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -16}, {40, -16}, {40, -16}}, color = {0, 0, 127}));
    connect(importWeatherModified.port_a, wall_North.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -12}, {40, -12}, {40, -10}}, color = {191, 0, 0}));
    connect(importWeatherModified.tworadiationport1, wall_East.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -48}, {40, -48}, {40, -48}}));
    connect(importWeatherModified.tworadiationport1, wall_North.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -8}, {40, -8}, {40, -8}}));
    connect(importWeatherModified.port_a, wall_East.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -52}, {40, -52}, {40, -50}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, roof.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, 24}, {40, 24}, {40, 24}}, color = {0, 0, 127}));
    connect(importWeatherModified.tworadiationport1, roof.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, 32}, {40, 32}, {40, 32}}));
    connect(importWeatherModified.port_a, roof.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, 28}, {40, 28}, {40, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_South.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, 30}, {-60, 30}, {-60, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_west.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -10}, {-60, -10}, {-60, -10}}, color = {191, 0, 0}));
    connect(radExchange.north, wall_North.star1) annotation(
      Line(points = {{-14, 20}, {-14, 36}, {6, 36}, {6, -4}, {20, -4}}, color = {191, 0, 0}));
    connect(radExchange.west, wall_west.star1) annotation(
      Line(points = {{-20, 10}, {-28, 10}, {-28, -4}, {-40, -4}, {-40, -4}}, color = {191, 0, 0}));
    connect(radExchange.floor, floor.star1) annotation(
      Line(points = {{-4, 0}, {-4, 0}, {-4, -10}, {-28, -10}, {-28, -44}, {-40, -44}, {-40, -44}}, color = {191, 0, 0}));
    connect(radExchange.ceiling, roof.star1) annotation(
      Line(points = {{-4, 20}, {6, 20}, {6, 36}, {20, 36}, {20, 36}, {20, 36}}, color = {191, 0, 0}));
    connect(radExchange.east, wall_East.star1) annotation(
      Line(points = {{0, 10}, {6, 10}, {6, -44}, {20, -44}, {20, -44}}, color = {191, 0, 0}));
    connect(radExchange.south, wall_South.star1) annotation(
      Line(points = {{-14, 0}, {-14, 0}, {-14, -4}, {-28, -4}, {-28, 36}, {-40, 36}, {-40, 36}}, color = {191, 0, 0}));
    connect(window_Status_hour_of_day_Winter.HourOfDay, importWeatherModified.HourOfDay) annotation(
      Line(points = {{38, 74}, {56, 74}, {56, 56}, {76, 56}, {76, 60}, {76, 60}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Winter.Window_Status, wall_North.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -2}, {40, -2}, {40, 0}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Winter.Window_Status, wall_East.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -40}, {40, -40}, {40, -40}}, color = {255, 127, 0}));
  protected
    annotation(
      experiment(StartTime = 0, StopTime = 5.0976e+06, Tolerance = 0.001, Interval = 10));
  end PrototypeSampleCase_SaintTropez_WinterJanFeb;

  model PrototypeSampleCase_SaintTropez_WinterDec
    HomeAssignment_2023.Wall wall_North(A_window = 2.25, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = 180, h = 3.2, horizontal = false, l = 4, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_East(A_window = 9, K = 3.5, N_pane = 1.526, U = 1.5, alpha_SW = 0.95, beta = 90, c = {921, 850, 850, 1100}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_exterior, d = {0.10, 0.25, 0.20, 0.01}, d_pane = 4, eps_LW = 0.95, frameratio = 0.3, gamma = -90, h = 3.2, horizontal = false, l = 6, lambda = {1.333, 0.043, 0.727, 0.415}, lat = 43.273, lowestFloor = false, n = 4, number_panes = 3, outside = true, rho(each displayUnit = "kg/m3") = {2000, 50, 1900, 1248}, rho_floor = 0.2, roof = false, withWindow = true) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_South(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 0, h = 3.2, horizontal = false, l = 4, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall wall_west(alpha_SW = 0, beta = 90, c = {1300, 850, 1300}, control_type_down_or_in = HeatConvUniversal.vert_interior, control_type_up_or_out = HeatConvUniversal.vert_interior, d = {0.01, 0.12, 0.01}, eps_LW = 0.95, gamma = 90, h = 3.2, horizontal = false, l = 6, lambda = {1.4, 0.04, 1.4}, lat = 43.273, lowestFloor = false, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {750, 60, 750}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall roof(alpha_SW = 0.95, beta = 0, c = {950, 900, 1100}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.25, 0.10, 0.01}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {2.0, 0.043, 0.415}, lat = 43.273, lowestFloor = false, n = 3, outside = true, rho(each displayUnit = "kg/m3") = {2400, 100, 1248}, rho_floor = 0.2, roof = true, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {30, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    HomeAssignment_2023.Wall floor(T0(displayUnit = "K") = 289.15, alpha_SW = 0, beta = 0, c = {1600, 900, 950}, control_type_down_or_in = HeatConvUniversal.hor_up, control_type_up_or_out = HeatConvUniversal.hor_down, d = {0.021, 0.12, 0.30}, eps_LW = 0.95, gamma = 0, h = 6, horizontal = true, l = 4, lambda = {0.13, 0.043, 2.0}, lat = 43.273, lowestFloor = true, n = 3, outside = false, rho(each displayUnit = "kg/m3") = {462, 100, 2400}, rho_floor = 0.2, roof = false, withWindow = false) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.SolarPositionModel solarPositionModel(L = 43.273) annotation(
      Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow = 0.01, T_ref = 289.15, alpha = 1) annotation(
      Placement(visible = true, transformation(origin = {-90, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    BaseLib.Airload airload(T(fixed = true), T0(displayUnit = "K") = 289.15, V = 76.8, c = 1007, rho = 1.19) annotation(
      Placement(visible = true, transformation(origin = {-10, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(visible = true, transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-12, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ImportWeatherModified importWeatherModified annotation(
      Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature fixedTemperature1(T = 16) annotation(
      Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.RadExchange radExchange(epsC = 0.95, epsEast = 0.95, epsF = 0.95, epsNorth = 0.95, epsSouth = 0.95, epsWest = 0.95, height = 3.2, length = 6, width = 4) annotation(
      Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    HomeAssignment_2023.Window_Status_hour_of_day_Winter window_Status_hour_of_day_Winter annotation(
      Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(fixedHeatFlow.port, floor.port_a) annotation(
      Line(points = {{-80, -50}, {-60, -50}, {-60, -50}, {-60, -50}}, color = {191, 0, 0}));
    connect(solarPositionModel.solarposition, roof.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, 36}, {40, 36}, {40, 36}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_North.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -4}, {40, -4}, {40, -4}}, color = {0, 0, 127}));
    connect(solarPositionModel.solarposition, wall_East.solarposition) annotation(
      Line(points = {{90, 60}, {90, 60}, {90, -44}, {40, -44}, {40, -44}}, color = {0, 0, 127}));
    connect(wall_South.port_b, port_a) annotation(
      Line(points = {{-40, 30}, {-34, 30}, {-34, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, floor.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -50}, {-40, -50}}, color = {191, 0, 0}));
    connect(airload.port_a, port_a) annotation(
      Line(points = {{-10, -60}, {-10, -60}, {-10, -30}, {-10, -30}}, color = {191, 0, 0}));
    connect(port_a, wall_west.port_b) annotation(
      Line(points = {{-10, -30}, {-34, -30}, {-34, -10}, {-40, -10}}, color = {191, 0, 0}));
    connect(port_a, wall_East.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -50}, {20, -50}, {20, -50}, {20, -50}}, color = {191, 0, 0}));
    connect(port_a, wall_North.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, -10}, {20, -10}, {20, -10}}, color = {191, 0, 0}));
    connect(port_a, roof.port_b) annotation(
      Line(points = {{-10, -30}, {12, -30}, {12, 30}, {20, 30}, {20, 30}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, wall_East.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -56}, {40, -56}, {40, -56}}, color = {0, 0, 127}));
    connect(importWeatherModified.windspeed, wall_North.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, -16}, {40, -16}, {40, -16}}, color = {0, 0, 127}));
    connect(importWeatherModified.port_a, wall_North.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -12}, {40, -12}, {40, -10}}, color = {191, 0, 0}));
    connect(importWeatherModified.tworadiationport1, wall_East.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -48}, {40, -48}, {40, -48}}));
    connect(importWeatherModified.tworadiationport1, wall_North.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, -8}, {40, -8}, {40, -8}}));
    connect(importWeatherModified.port_a, wall_East.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, -52}, {40, -52}, {40, -50}}, color = {191, 0, 0}));
    connect(importWeatherModified.windspeed, roof.u) annotation(
      Line(points = {{74, 60}, {72, 60}, {72, 24}, {40, 24}, {40, 24}}, color = {0, 0, 127}));
    connect(importWeatherModified.tworadiationport1, roof.tworadiationport1) annotation(
      Line(points = {{68, 60}, {68, 60}, {68, 32}, {40, 32}, {40, 32}}));
    connect(importWeatherModified.port_a, roof.port_a) annotation(
      Line(points = {{64, 60}, {64, 60}, {64, 28}, {40, 28}, {40, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_South.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, 30}, {-60, 30}, {-60, 30}}, color = {191, 0, 0}));
    connect(fixedTemperature1.port, wall_west.port_a) annotation(
      Line(points = {{-80, 10}, {-70, 10}, {-70, -10}, {-60, -10}, {-60, -10}}, color = {191, 0, 0}));
    connect(radExchange.north, wall_North.star1) annotation(
      Line(points = {{-14, 20}, {-14, 36}, {6, 36}, {6, -4}, {20, -4}}, color = {191, 0, 0}));
    connect(radExchange.west, wall_west.star1) annotation(
      Line(points = {{-20, 10}, {-28, 10}, {-28, -4}, {-40, -4}, {-40, -4}}, color = {191, 0, 0}));
    connect(radExchange.floor, floor.star1) annotation(
      Line(points = {{-4, 0}, {-4, 0}, {-4, -10}, {-28, -10}, {-28, -44}, {-40, -44}, {-40, -44}}, color = {191, 0, 0}));
    connect(radExchange.ceiling, roof.star1) annotation(
      Line(points = {{-4, 20}, {6, 20}, {6, 36}, {20, 36}, {20, 36}, {20, 36}}, color = {191, 0, 0}));
    connect(radExchange.east, wall_East.star1) annotation(
      Line(points = {{0, 10}, {6, 10}, {6, -44}, {20, -44}, {20, -44}}, color = {191, 0, 0}));
    connect(radExchange.south, wall_South.star1) annotation(
      Line(points = {{-14, 0}, {-14, 0}, {-14, -4}, {-28, -4}, {-28, 36}, {-40, 36}, {-40, 36}}, color = {191, 0, 0}));
    connect(window_Status_hour_of_day_Winter.HourOfDay, importWeatherModified.HourOfDay) annotation(
      Line(points = {{38, 74}, {56, 74}, {56, 56}, {76, 56}, {76, 60}, {76, 60}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Winter.Window_Status, wall_North.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -2}, {40, -2}, {40, 0}}, color = {255, 127, 0}));
    connect(window_Status_hour_of_day_Winter.Window_Status, wall_East.WindowOpeningInput) annotation(
      Line(points = {{40, 70}, {52, 70}, {52, -40}, {40, -40}, {40, -40}}, color = {255, 127, 0}));
  protected
    annotation(
      experiment(StartTime = 2.88576e+07, StopTime = 3.1536e+07, Tolerance = 0.001, Interval = 10));
  end PrototypeSampleCase_SaintTropez_WinterDec;
  annotation(
    uses(Modelica(version = "3.2.3")));
end HomeAssignment_2023;
