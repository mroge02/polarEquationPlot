# polarEquationPlot
Mathematica package for plotting a level set given as an equation in polar coordinates

This is a low priority for me at present. Someday if I have the time, I'll make a full paclet. Feel free to let me know of problems. When I have time, I will try to fix them.

**Examples**

    polarEquationPlot[r^2 == t, {t, 0, 60}, {r, -8, 8}]

    polarEquationPlot[(r - 2)^2 + t^2 == 9, {t, -Pi, Pi}, {r, -1, 5}, 
     PolarAxes -> True, PolarGridLines -> Automatic]

    GraphicsColumn[{
      polarEquationPlot[r^2 == t, {t, 0, 60}, {r, -8, 8}, 
       PolarAxes -> Automatic, Axes -> True, PolarGridLines -> Automatic, 
       ImageSize -> Large],
      getLastPolarMeshPlot[] (* draws underlying mesh of last plot *)
      }]

    polarEquationPlot[
     10 Cos[(r - 1)]^2 + t^3 == 9, {t, -2.5, 2.5}, {r, -6, 6}, 
     MaxCellMeasure -> {"Curvature" -> 0.05}, PolarAxes -> Automatic, 
     Axes -> True, PolarGridLines -> Automatic, 
     GridLinesStyle -> RGBColor[0.85, 0.8, 0.75]]

    lastPolarEquationPlotData (* data used to generate last plot *)
