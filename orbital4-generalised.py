import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# physics engine

def force_gravity(pos1, pos2, m1, m2, G=1.0):
    """
    Calculate gravitational force on body 1 due to body 2.
    
    Inputs:
        pos1, pos2: 3D position vectors [x, y, z]
        m1, m2: masses
        G: gravitational constant (default 1.0 for dimensionless units)
    
    Returns:
        force_vector: 3D force on body 1
    """
    r_vec = pos2 - pos1
    r_mag = np.linalg.norm(r_vec)
    
    if r_mag == 0:
        return np.zeros(3)
    
    # F = G * m1 * m2 / r^2, direction: towards body 2
    force_mag = G * m1 * m2 / (r_mag ** 2)
    force_vec = force_mag * (r_vec / r_mag)
    
    return force_vec


def calculate_accelerations(bodies, G=1.0):
    """
    Calculate accelerations for all bodies due to mutual gravitational attraction.
    
    Inputs:
        bodies: list of dicts with keys "mass", "pos", "vel"
        G: gravitational constant
    
    Returns:
        accelerations: list of 3D acceleration vectors
    """
    n = len(bodies)
    accelerations = [np.zeros(3) for _ in range(n)]
    
    # O(N^2) force calculation - all pairs
    for i in range(n):
        for j in range(n):
            if i != j:
                force = force_gravity(
                    bodies[i]["pos"], 
                    bodies[j]["pos"], 
                    bodies[i]["mass"], 
                    bodies[j]["mass"], 
                    G
                )
                # a = F / m
                accelerations[i] += force / bodies[i]["mass"]
    
    return accelerations


def integrate_euler(bodies, dt, G=1.0):
    """
    Euler integration: update positions and velocities for one timestep.
    
    Inputs:
        bodies: list of body dicts (modified in place)
        dt: timestep
        G: gravitational constant
    """
    # Calculate accelerations at current state
    accelerations = calculate_accelerations(bodies, G)
    
    # Update velocities and positions
    for i, body in enumerate(bodies):
        body["vel"] = body["vel"] + accelerations[i] * dt
        body["pos"] = body["pos"] + body["vel"] * dt


def integrate_leapfrog(bodies, dt, G=1.0):
    """
    Leapfrog (Verlet) integration: symplectic integrator with better energy conservation.
    
    Second-order method, much better for orbital mechanics than Euler.
    
    Algorithm (Kick-Drift-Kick):
    1. Half-(time)step velocity update
    2. Full-step position update  
    3. Half-step velocity update
    
    Inputs:
        bodies: list of body dicts (modified in place)
        dt: timestep
        G: gravitational constant
    """
    # First half-kick:
    accelerations = calculate_accelerations(bodies, G)
    for i, body in enumerate(bodies):
        body["vel"] = body["vel"] + accelerations[i] * (dt / 2.0)
    
    # Drift:
    for body in bodies:
        body["pos"] = body["pos"] + body["vel"] * dt
    
    # Second half-kick:
    accelerations = calculate_accelerations(bodies, G)
    for i, body in enumerate(bodies):
        body["vel"] = body["vel"] + accelerations[i] * (dt / 2.0)


def calculate_timestep(bodies, target_steps_per_orbit=100):
    """
    Auto-calculate appropriate timestep based on initial conditions.
    
    Estimates smallest orbital period and divides by target steps.
    
    Inputs:
        bodies: list of body dicts
        target_steps_per_orbit: desired number of timesteps per orbit
    
    Returns:
        dt: timestep
    """
    min_period = float("inf")
    
    # Estimate orbital periods for each pair
    for i in range(len(bodies)):
        for j in range(i + 1, len(bodies)):
            r_vec = bodies[j]["pos"] - bodies[i]["pos"]
            r = np.linalg.norm(r_vec)
            
            if r == 0:
                continue
            
            # For circular orbit: v = sqrt(G*M/r), period = 2*pi*r/v
            # use total mass for two-body approximation
            M_total = bodies[i]["mass"] + bodies[j]["mass"]
            v_circular = np.sqrt(M_total / r)  # G=1
            period = 2 * np.pi * r / v_circular
            
            min_period = min(min_period, period)
    
    if min_period == float("inf"):
        # Fallback if no valid period found
        return 0.01
    
    dt = min_period / target_steps_per_orbit
    return dt


# simulation state

def initialise_two_body(eccentricity=0.0):
    """
    Initialise a 2-body system: one massive central body, one orbiting body.
    
    Args:
        eccentricity: Orbital eccentricity (0 = circular, 0 < e < 1 = elliptical)
    
    Returns:
        bodies: list of body dicts
    """
    # Body 1: Central, massive body
    # For simplicity, keep it at origin with zero velocity
    body1 = {
        "mass": 1.0,
        "pos": np.array([0.0, 0.0, 0.0]),
        "vel": np.array([0.0, 0.0, 0.0]),
        "trail": []
    }
    
    # Body 2: Orbiting body
    # Place at periapsis (closest approach)
    periapsis = 10.0
    M_central = body1["mass"]
    
    # Velocity at periapsis: v² = GM*(1+e)/r_peri
    v_periapsis = np.sqrt(M_central * (1 + eccentricity) / periapsis)
    
    body2 = {
        "mass": 0.01,
        "pos": np.array([periapsis, 0.0, 0.0]),
        "vel": np.array([0.0, v_periapsis, 0.0]),
        "trail": []
    }
    
    return [body1, body2]


def initialise_binary_star(eccentricity=0.0):
    """
    Initialise a binary star system with unequal masses orbiting their barycenter.
    
    Args:
        eccentricity: Orbital eccentricity (0 = circular, 0 < e < 1 = elliptical)
    
    Returns:
        bodies: list of body dicts
    """
    # Binary system with unequal masses for visually distinct orbits
    mass1 = 2.0  # More massive star
    mass2 = 1.0  # Less massive star
    M_total = mass1 + mass2
    
    # Define the orbit at periapsis (closest approach)
    periapsis_separation = 4.0  # Closest distance between stars
    
    # Semi-major axis of relative orbit: r_peri = a(1-e)
    # For circular orbit (e=0): a = r_peri
    a_relative = periapsis_separation / (1 - eccentricity) if eccentricity < 1.0 else periapsis_separation
    
    # Distances from COM at periapsis (preserve COM at origin)
    r1 = periapsis_separation * mass2 / M_total
    r2 = periapsis_separation * mass1 / M_total
    
    # Velocity at periapsis using vis-viva equation
    # For circular (e=0): v² = GM/r
    # For elliptical: v² = GM*(1+e)/(a(1-e)) = GM*(1+e)/r_peri
    v_relative = np.sqrt(M_total * (1 + eccentricity) / periapsis_separation)
    
    # Individual velocities (conserving momentum)
    v1 = v_relative * mass2 / M_total
    v2 = v_relative * mass1 / M_total
    
    star1 = {
        "mass": mass1,
        "pos": np.array([-r1, 0.0, 0.0]),
        "vel": np.array([0.0, -v1, 0.0]),
        "trail": []
    }
    
    star2 = {
        "mass": mass2,
        "pos": np.array([r2, 0.0, 0.0]),
        "vel": np.array([0.0, v2, 0.0]),
        "trail": []
    }
    
    return [star1, star2]


def initialise_sun_earth_moon(earth_eccentricity=0.0, moon_eccentricity=0.0):
    """
    Initialise Sun-Earth-Moon system (scaled dimensionless units).
    
    Args:
        earth_eccentricity: Eccentricity of Earth's orbit around Sun
        moon_eccentricity: Eccentricity of Moon's orbit around Earth
    
    Mass ratios approximately:
    - Sun: 1.0
    - Earth: 0.000003 (actually ~3e-6, using 0.003 for visibility)
    - Moon: 0.0000001 (actually ~3.7e-8, using 0.0001 for visibility)
    
    Returns:
        bodies: list of body dicts
    """
    # Sun at origin
    sun = {
        "mass": 1.0,
        "pos": np.array([0.0, 0.0, 0.0]),
        "vel": np.array([0.0, 0.0, 0.0]),
        "trail": []
    }
    
    # Earth orbiting Sun (at periapsis)
    earth_periapsis = 10.0
    v_earth = np.sqrt(sun["mass"] * (1 + earth_eccentricity) / earth_periapsis)
    
    earth = {
        "mass": 0.003,  # Scaled up for visible wobble
        "pos": np.array([earth_periapsis, 0.0, 0.0]),
        "vel": np.array([0.0, v_earth, 0.0]),
        "trail": []
    }
    
    # Moon orbiting Earth (at periapsis, relative to Earth's position/velocity)
    moon_periapsis = 0.3  # Distance from Earth
    # Moon's velocity relative to Earth for orbit around Earth
    v_moon_rel = np.sqrt(earth["mass"] * (1 + moon_eccentricity) / moon_periapsis)
    
    moon = {
        "mass": 0.0001,
        "pos": earth["pos"] + np.array([moon_periapsis, 0.0, 0.0]),
        "vel": earth["vel"] + np.array([0.0, v_moon_rel, 0.0]),
        "trail": []
    }
    
    return [sun, earth, moon]


def initialise_figure_eight():
    """
    Initialise the famous figure-eight three-body solution.
    
    This is a special periodic solution where three equal masses
    chase each other in a figure-eight pattern.
    
    Returns:
        bodies: list of body dicts
    """
    mass = 1.0
    
    # Known initial conditions for figure-eight orbit
    # (Chenciner & Montgomery, 2000)
    body1 = {
        "mass": mass,
        "pos": np.array([-0.97000436, 0.24308753, 0.0]),
        "vel": np.array([0.4662036850, 0.4323657300, 0.0]),
        "trail": []
    }
    
    body2 = {
        "mass": mass,
        "pos": np.array([0.0, 0.0, 0.0]),
        "vel": np.array([-0.93240737, -0.86473146, 0.0]),
        "trail": []
    }
    
    body3 = {
        "mass": mass,
        "pos": np.array([0.97000436, -0.24308753, 0.0]),
        "vel": np.array([0.4662036850, 0.4323657300, 0.0]),
        "trail": []
    }
    
    return [body1, body2, body3]


def initialise_hierarchical_triple(binary_eccentricity=0.0, third_body_eccentricity=0.0):
    """
    Initialise a hierarchical triple system: binary pair + distant third body.
    
    Args:
        binary_eccentricity: Eccentricity of the close binary pair
        third_body_eccentricity: Eccentricity of the third body's orbit around the binary
    
    Returns:
        bodies: list of body dicts
    """
    # Close binary pair with equal masses
    mass_binary = 1.0
    M_binary_total = 2 * mass_binary
    
    # Binary at periapsis
    binary_periapsis = 2.0
    r_binary = binary_periapsis / 2.0  # Each star's distance from binary COM
    
    # Binary velocity at periapsis
    v_rel_binary = np.sqrt(M_binary_total * (1 + binary_eccentricity) / binary_periapsis)
    v_binary = v_rel_binary / 2.0  # Each star gets half the relative velocity
    
    body1 = {
        "mass": mass_binary,
        "pos": np.array([-r_binary, 0.0, 0.0]),
        "vel": np.array([0.0, -v_binary, 0.0]),
        "trail": []
    }
    
    body2 = {
        "mass": mass_binary,
        "pos": np.array([r_binary, 0.0, 0.0]),
        "vel": np.array([0.0, v_binary, 0.0]),
        "trail": []
    }
    
    # Distant third body orbiting the binary's barycenter
    mass_third = 0.5
    third_periapsis = 15.0
    
    # Third body velocity at periapsis
    v_third = np.sqrt(M_binary_total * (1 + third_body_eccentricity) / third_periapsis)
    
    # CRITICAL: Conserve momentum! Binary COM must move to balance third body
    v_correction = -mass_third * np.array([-v_third, 0.0, 0.0]) / M_binary_total
    body1["vel"] = body1["vel"] + v_correction
    body2["vel"] = body2["vel"] + v_correction
    
    body3 = {
        "mass": mass_third,
        "pos": np.array([0.0, third_periapsis, 0.0]),
        "vel": np.array([-v_third, 0.0, 0.0]),
        "trail": []
    }
    
    return [body1, body2, body3]


# visualisation

def setup_plot():
    """
    Set up 3D matplotlib figure and axis.
    
    Returns:
        fig, ax: matplotlib figure and 3D axis
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    
    # Set axis limits
    limit = 15
    ax.set_xlim([-limit, limit])
    ax.set_ylim([-limit, limit])
    ax.set_zlim([-limit, limit])
    
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("N-Body Gravitational Simulation")
    
    return fig, ax


def update_trails(bodies, max_trail_length=5000):
    """
    Update trail buffers for each body (circular buffer).
    
    Args:
        bodies: list of body dicts
        max_trail_length: maximum number of points to keep
    """
    for body in bodies:
        body["trail"].append(body["pos"].copy())
        
        # Keep only last max_trail_length points
        if len(body["trail"]) > max_trail_length:
            body["trail"].pop(0)


# animation

class NBodyAnimation:
    def __init__(self, bodies, dt, steps_per_frame=10, reference_body_idx=None, integrator="leapfrog"):
        """
        Initialise animation.
        
        Args:
            bodies: list of body dicts
            dt: timestep
            steps_per_frame: physics steps between visualization updates
            reference_body_idx: index of body to use as reference frame (None for global)
            integrator: 'euler' or 'leapfrog' (default: 'leapfrog')
        """
        self.bodies = bodies
        self.dt = dt
        self.steps_per_frame = steps_per_frame
        self.reference_body_idx = reference_body_idx
        self.frame_count = 0
        
        # Select integration method
        if integrator == "leapfrog":
            self.integrate = integrate_leapfrog
        else:
            self.integrate = integrate_euler
        
        print(f"Using {integrator} integration")
        
        # Set up plot
        self.fig, self.ax = setup_plot()
        
        # Initialise plot elements
        self.body_plots = []
        self.trail_plots = []
        
        colors = ["gold", "blue", "red", "green", "purple"]
        
        for i, body in enumerate(self.bodies):
            # Body marker
            body_plot, = self.ax.plot(
                [body["pos"][0]], 
                [body["pos"][1]], 
                [body["pos"][2]], 
                "o", 
                color=colors[i % len(colors)],
                markersize=10,
                label=f"Body {i+1}"
            )
            self.body_plots.append(body_plot)
            
            # Trail line
            trail_plot, = self.ax.plot(
                [], [], [], 
                '-', 
                color=colors[i % len(colors)],
                alpha=0.6,
                linewidth=1
            )
            self.trail_plots.append(trail_plot)
        
        self.ax.legend()
    
    def update(self, frame):
        """
        Animation update function called each frame.
        """
        # Perform multiple physics steps per frame
        for _ in range(self.steps_per_frame):
            self.integrate(self.bodies, self.dt)
            update_trails(self.bodies)
        
        # Get reference position if reference frame is set
        ref_pos = np.zeros(3)
        if self.reference_body_idx is not None:
            ref_pos = self.bodies[self.reference_body_idx]["pos"]
        
        # Update visualization
        for i, body in enumerate(self.bodies):
            # Transform to reference frame
            display_pos = body["pos"] - ref_pos
            
            # Update body position
            self.body_plots[i].set_data([display_pos[0]], [display_pos[1]])
            self.body_plots[i].set_3d_properties([display_pos[2]])
            
            # Update trail
            if len(body["trail"]) > 0:
                trail_array = np.array(body["trail"])
                # Transform trail to reference frame
                display_trail = trail_array - ref_pos
                self.trail_plots[i].set_data(display_trail[:, 0], display_trail[:, 1])
                self.trail_plots[i].set_3d_properties(display_trail[:, 2])
        
        self.frame_count += 1
        if self.frame_count % 100 == 0:
            print(f"Frame {self.frame_count}, Time: {self.frame_count * self.steps_per_frame * self.dt:.2f}")
        
        return self.body_plots + self.trail_plots
    
    def run(self):
        """
        Start the animation.
        """
        anim = FuncAnimation(
            self.fig, 
            self.update, 
            interval=20,  # ms between frames (~50 FPS)
            blit=False,   # blit=True doesn't work well with 3D
            cache_frame_data=False
        )
        plt.show()


# main

if __name__ == "__main__":    
    # Choose initialization function for different scenarios:
    # - initialise_two_body(e)                 : Simple star-planet system
    # - initialise_binary_star(e)              : Unequal mass binary stars
    # - initialise_sun_earth_moon(e1, e2)      : 3-body hierarchical system
    # - initialise_figure_eight()              : Famous figure-8 periodic orbit (no eccentricity param)
    # - initialise_hierarchical_triple(e1, e2) : Binary pair + distant third body

    # Eccentricity (e):
    #  e = 0.0 → Perfect circle
    #  e = 0.3 → Slightly elliptical
    #  e = 0.5 → Moderately elliptical
    #  e = 0.7 → Very elliptical
    #  e → 1.0 → Nearly parabolic (unbounded)
    
    # Examples:
    bodies = initialise_binary_star(eccentricity=0.5) # Elliptical binary
    # bodies = initialise_binary_star(eccentricity=0.0) # Circular binary
    # bodies = initialise_two_body(eccentricity=0.3) # Slightly elliptical planet orbit
    # bodies = initialise_sun_earth_moon(earth_eccentricity=0.0, moon_eccentricity=0.1)
    # bodies = initialise_hierarchical_triple(binary_eccentricity=0.2, third_body_eccentricity=0.4)
    # bodies = initialise_figure_eight() # No eccentricity control for special solutions
    # bodies = initialise_sun_earth_moon(earth_eccentricity=0.017, moon_eccentricity=0.055) # Sun-Earth-Moon system viewed in a reference frame
    
    # Reference frame:
    # - None: Global inertial frame (COM stays at origin)
    # - 0: Body 1's frame (Body 1 appears stationary)
    # - 1: Body 2's frame (Body 2 appears stationary)
    # - etc.
    reference_body_idx = None
    
    # Integration method: "euler" or "leapfrog"
    # Leapfrog is more accurate and conserves energy better for orbits
    integrator = "leapfrog"
    
    # Animation params:
    # - target_steps_per_orbit: More steps = better accuracy but slower
    # - steps_per_frame: More steps = faster simulation, fewer = smoother animation
    target_steps_per_orbit = 200 # Increased for better accuracy with Leapfrog
    steps_per_frame = 5
    
    # ========================================================================
    
    print(f"Initialising {len(bodies)}-body system")
    print(f"Reference frame: {"Global" if reference_body_idx is None else f"Body {reference_body_idx + 1}"}")
    
    # Calculate timestep
    dt = calculate_timestep(bodies, target_steps_per_orbit=target_steps_per_orbit)
    print(f"Calculated timestep: {dt:.6f}")
        
    animation = NBodyAnimation(bodies, dt, steps_per_frame=steps_per_frame, reference_body_idx=reference_body_idx, integrator=integrator)
    animation.run()
