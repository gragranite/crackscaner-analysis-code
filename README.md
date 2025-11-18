Crackscaner Ver.1.0 –README 

© Takato TAKEMURA, Geomechnics Lab., Nihon University, Tokyo, Japan

1. Overview
This document describes the fully updated functionality of the crack-scan analyzer:
svg_viewer_boundary_scale_axes_choice_axisflip_scanline_multi_angle.py

The tool provides:
• SVG loading and geometric parsing
• Automatic extraction of crack segments and closed boundaries
• Pixel-to-mm scaling and unit conversion
• Coordinate-axis flipping (left–right, up–down)
• Multi-tab Tkinter GUI
• Single and multiple horizontal scanline analysis
• Full 0–180° angle sweep with multiple scanlines per angle
• Tab-based display for each angle
• Convergence plots (separate window)
• Rose diagram (0–360°) with adjustable radius and CSV/SVG export
• Crack-length histogram (normalized probability)
• Full CSV export of all numerical results
2. Input Data Requirements
(1) SVG format
Files should be exported from Adobe Illustrator or equivalent, containing line-based geometry.

(2) Closed boundary
A closed boundary (polyline, polygon, rect, or closed path) defines the analysis domain.
The boundary must be visually distinguishable (commonly red).

(3) Crack traces
Cracks must be drawn as line/polyline/path elements.
Curved paths are internally sampled.

(4) Scale
Users specify pixel-to-mm conversion in the control panel.
3. Control Panel and Functions
The GUI consists of regions:

A. File loading
[Open file...] loads SVG and displays geometry.

B. Scale & axes
- Unit → mm conversion
- Projection axis selection (XY/XZ/YZ)
- Scale bar length selection
- Axis flipping (horizontal/vertical)
- Mini-axis glyph auto-adjusting

C. Single scanline mode
Input: y-value or auto-center
Outputs:
• scanline length (units & mm)
• intersection count
• density (per unit, per mm)
• intersection markers

D. Multiple horizontal scanlines
Modes: Even spacing / Incremental
User sets: number of lines
Outputs:
• multiple scanlines covering domain
• cumulative statistics
• convergence plot (separate window)

E. Angle sweep (0–180°)
Settings:
• Step: 10° or 15°
• Multiple scanlines per angle (same as horizontal mode)
Outputs:
• Tabs for each angle
• Angle-specific convergence plot
• Scanline geometry & hits
4. Rose Diagram Module
Settings:
• L_total (total scanline length) override
• Max radius input (blank = auto)
• 0–360° direction using mirrored 0–180° data

Features:
• Density expressed as count/mm
• θ = 0° oriented upward
• Clockwise positive rotation
• Axis labels consistent with projection & flip settings

Exports:
• SVG rose diagram
• CSV containing: angle_deg, density_per_mm, total_length_used, mm_per_unit
5. Crack-Length Histogram
Function:
• Computes length of each SVG crack segment (units → mm)
• Produces normalized histogram (sum = 1)

Outputs:
• Histogram plot (Tk window)
• SVG histogram
• CSV containing: index, length_units, length_mm
6. CSV Export Specification
Full CSV export includes all modes (single, multiple, angle sweep).

Columns:
angle_deg,
mode,
direction,
index,
scanline_position,
length_units,
length_mm,
hit_count,
density_per_unit,
density_per_mm,
cumulative_length,
cumulative_hits,
cumulative_density

This enables statistical analysis in Python/R/Excel.
7. Execution Environment
• macOS 
• Python ≥ 3.x
Required libraries:
• tkinter
• numpy
• matplotlib
• standard library (no third-party dependencies)
8. Notes and Limitations
• Curved cracks are approximated via path sampling.
• Boundary must be fully closed or no scanlines are valid.
• Extremely large SVGs may reduce performance.
• Illustrator paths should be constructed carefully for accuracy.
9. References
Oda, M., Katsube, T., Takemura, T. (2002). Microcrack evolution and brittle failure of Inada granite in triaxial compression tests at 140 MPa. Journal of Geophysical Research, 107(B10), 2233.
Takemura, T., Oda, M. (2004). Stereology-based fabric analysis of microcracks in damaged granite. Tectonophysics, 387, 131–150.
<img width="432" height="647" alt="image" src="https://github.com/user-attachments/assets/ad4e9770-5a8d-4acd-b9f5-807d75af06b8" />
