# N-Body Gravitational Simulation - Project Specification

## Project Overview
Interactive 3D simulation of gravitational n-body systems using matplotlib. Focus on 2-5 bodies with physically informative visualisations of orbital dynamics.

---

## Simulation Notes

### Target Systems to simulate
- Sun-Earth-Moon
- Binary star systems
- Star-planet systems
- Star with 2 planets
- Planet with 2 stars
- General: 2-5 bodies with arbitrary masses

### Physics Framework
- **Coordinate system**: Global inertial frame (no hardcoded center)
- **Forces**: Gravitational attraction between all body pairs
- **Units**: Dimensionless (masses, distances, velocities are relative)
- **Integration**: Start with Euler method, migrate to Leapfrog, keep option for RK4 later
- **Timestep**: Fixed timestep auto-calculated from initial conditions (Option A)
  - Calculate once at initialisation based on system properties
  - Target: ~100-200 steps per smallest expected orbital period
  - Use throughout entire simulation run

### Simulation Behavior
- **Duration**: Runs until user manually stops (Reset/Stop button) or collision detected
- **Collision detection**: Auto-stop when body centers overlap (distance ≤ sum of radii or small threshold)
- **Escape handling**: No auto-stop for escapes - let simulation continue, user observes and stops manually
- **Frame rate**: Update visualisation every 5-10 physics steps (~30-60 FPS effective)

---

## UI Design Notes

### Control Flow
1. **Set** button: Initialise system with current slider values
2. **Run/Pause** button: Start/pause animation (toggle)
3. **Reset/Stop** button: Stop simulation, return to initial conditions setup mode

**Key rule**: No parameter changes allowed during simulation (even when paused). Sliders only set initial conditions.

### Widget Layout
```
[Top Control Panel]
├─ Number of bodies: [dropdown: 2, 3, 4, 5]
├─ [Set] [Run/Pause] [Reset/Stop]

[Body Parameters - Grouped by Body]
Body 1:  
  Mass [log slider: 0.001 to 1000] {with markers}
  Position X [linear slider: -50 to 50]
  Position Y [linear slider: -50 to 50]
  Position Z [linear slider: -50 to 50]
  Velocity X [linear slider: -10 to 10]
  Velocity Y [linear slider: -10 to 10]
  Velocity Z [linear slider: -10 to 10]
         
Body 2:  
  [same slider structure]
  
[... repeat for bodies 3-5 as needed]
```

**Total sliders**: 7 per body × N bodies (14 for 2 bodies, 35 for 5 bodies)

### Slider Markers (Reference Values)
- **Purpose**: Help user set reasonable initial conditions
- **Update timing**: Only when sliders are being adjusted (not real-time during simulation)
- **Types of markers**:
  - **Escape velocity** (always shown): Labeled with originating body name (e.g., "Body_1")
  - **Hill sphere radius** (3+ bodies only): Labeled with body name
- **Multiple markers**: If multiple bodies affect a parameter, show all relevant markers with labels
- **Marker disappearance**: Hill sphere markers disappear if body count < 3

### Slider Scale Types
- **Mass**: Logarithmic scale [0.001, 1000] - covers 6 orders of magnitude
- **Position (X, Y, Z)**: Linear scale [-50, 50]
- **Velocity (X, Y, Z)**: Linear scale [-10, 10]

---

## Visualisation

### 3D Display
- **Library**: Matplotlib (mplot3d, animation)
- **View**: 3D plot with rotation capability
- **Bodies**: Rendered as spheres/points with labels
- **Trails**: Orbital paths shown behind each body

### Trail Rendering
- **Type**: Disappearing trails (fixed-size circular buffer)
- **Resolution**: Adaptively calculated by program
  - Target: ~200-500 trail points per orbit for helical detail
  - Falls back to "every Nth frame" if orbital period unknown
- **Buffer size**: Last 5000 points per body (prevents unbounded memory growth)
- **Rationale**: Must capture features like moon's helical path around planet orbiting star

### Performance Considerations
- Force calculations: O(N²) acceptable for N ≤ 5
- Trail rendering: Circular buffer, downsampled display
- Physics steps: Not every step needs visualisation update
- Matplotlib 3D limitations acknowledged, optimise within constraints

---

## Initial Conditions (Reset State)

When app starts or user clicks Reset:

**For 2 bodies:**
- Body 1: Mass = 1.0, Position = (0, 0, 0), Velocity = (0, 0, 0)
- Body 2: Mass = 0.01, Position = (10, 0, 0), Velocity = (0, v_circular, 0)
  - Where v_circular is calculated for stable circular orbit

**For 3+ bodies:**
- Arrange in "Copernican" configuration (lined up along x-axis)
- Body 1 (central/largest): Position = (0, 0, 0), Velocity = (0, 0, 0)
- Body 2, 3, ...: Position = (10, 0, 0), (20, 0, 0), (30, 0, 0), ...
- Each outer body given circular orbit velocity around Body 1
- Mass ratios: decreasing (e.g., 1.0, 0.1, 0.01, ...)

---

## Implementation Roadmap

### Minimal Working Prototype
- 2 bodies (fixed count)
- Euler integration with auto-calculated fixed dt
- 3D visualisation with adaptive disappearing trails
- All sliders: Mass (log), Position XYZ, Velocity XYZ (14 total)
- Buttons: Set, Run/Pause, Reset/Stop
- Escape velocity markers on sliders
- Auto-stop on collision
- Default reset state: stable 2-body circular orbit

### Core Features
- Variable body count (2-5) with dropdown
- Hill sphere markers (3+ bodies)
- Switch to Leapfrog integration (keep Euler as option)
- Improved default configurations for 3-5 bodies
- Better visualisation (body sizes proportional to mass, colors, labels)

### Possible Advanced Features
- RK4 integration option
- Save/load initial conditions to file
- Export animation capabilities
- Collision visualisation (highlight colliding bodies)
- Console output for system energy, momentum
- Performance profiling and optimisation
- Better UI aesthetics
- Optional: Qt/tkinter control panel

### Design Principles
- **Modularity**: Physics functions independent of UI/visualisation
- **Reusability**: Integration method swappable by changing one function call
- **Widget flexibility**: Keep widget code modular for potential Qt/tkinter migration
- **Pure functions**: Physics/reference calculations are stateless where possible

---

## Technical Constraints & Notes

### Python Environment
- Platform: macOS
- Python version: 3.13.11
- Environment: Conda
- Primary simulation library: Matplotlib (pyplot, animation, widgets, mplot3d)

### Limitations
- Matplotlib 3D not fastest for complex animations
- Widget aesthetics limited compared to Qt/tkinter
- Large N (>5) not prioritised
- No relativistic effects, no tidal forces

### Physical Accuracy Goals
- Moderate accuracy: informative and representative
- Not research-grade, but better than "high school level"
- Energy conservation depends on integration method quality
- Acceptable for qualitative orbit analysis and exploration

## Future Considerations
- Duration control slider
- Real-time parameter adjustment
- Infinite trail option (non-disappearing)
- Tidal forces, Roche limits
- Energy/momentum conservation display
- Multiple integration methods selectable by user
- Export to video file
- Load real solar system data
- Relativistic corrections
- N > 5 optimisation