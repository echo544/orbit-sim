"""
Interactive 3D Three-Body Orbital Simulation
Standalone desktop application using mpl

- added "Hill Sphere" limit (blue line) to prevent sun from stripping the moon away from earth
- improved binary planet physics (m_1 + m_2) for velocity "hints"
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button
from mpl_toolkits.mplot3d import Axes3D

G = 1.0

class Body: # Celestial body with position, velocity, and mass
    def __init__(self, name, mass, position, velocity, color, size):
        self.name = name
        self.mass = mass
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.color = color
        self.size = size
        self.trail_x = [position[0]]
        self.trail_y = [position[1]]
        self.trail_z = [position[2]]
        
    def add_to_trail(self):
        self.trail_x.append(self.position[0])
        self.trail_y.append(self.position[1])
        self.trail_z.append(self.position[2])
        if len(self.trail_x) > 500:
            self.trail_x.pop(0)
            self.trail_y.pop(0)
            self.trail_z.pop(0)

def compute_force(body1, body2):
    r_vec = body2.position - body1.position
    r_mag = np.linalg.norm(r_vec)
    if r_mag < 0.01: return np.array([0.0, 0.0, 0.0])
    force_mag = G * body1.mass * body2.mass / (r_mag ** 2)
    return force_mag * r_vec / r_mag

def update_system(bodies, dt):
    n = len(bodies)
    accelerations = []
    for i in range(n):
        total_force = np.array([0.0, 0.0, 0.0])
        for j in range(n):
            if i != j:
                total_force += compute_force(bodies[i], bodies[j])
        accelerations.append(total_force / bodies[i].mass)
    
    for i in range(n):
        bodies[i].velocity += accelerations[i] * dt
        bodies[i].position += bodies[i].velocity * dt
        bodies[i].add_to_trail()

class OrbitSimulator:
    def __init__(self):
        self.fig = plt.figure(figsize=(14, 9))
        self.fig.canvas.manager.set_window_title("3-Body Orbital Simulator")
        
        self.ax = self.fig.add_subplot(121, projection="3d")
        self.ax.set_facecolor("black")
        self.fig.patch.set_facecolor("#1a1a1a")
        
        # Initial parameters
        self.sun_mass = 1000.0
        self.planet_mass = 2.0
        self.planet_distance = 2.0
        self.planet_velocity = 22.0
        self.moon_mass = 0.1
        self.moon_distance = 0.15 # reduced default to be inside Hill Sphere
        self.moon_velocity = 3.6
        self.dt = 0.001
        
        self.running = False
        self.bodies = []
        self.hint_lines = {}
        
        self.setup_controls()
        self.reset_simulation()
        
        self.anim = FuncAnimation(self.fig, self.animate, interval=20, blit=False)
        
    def setup_controls(self):
        panel_left = 0.55
        panel_width = 0.35
        slider_height = 0.02
        spacing = 0.035
        y_pos = 0.90
        
        # Sliders
        ax_sun = self.fig.add_axes([panel_left, y_pos, panel_width, slider_height])
        self.slider_sun = Slider(ax_sun, "Sun Mass", 100, 2000, valinit=self.sun_mass, 
                                 color="yellow", valstep=50, valfmt="%.0f")
        self.slider_sun.on_changed(self.update_sun_mass)
        
        y_pos -= spacing
        ax_planet_mass = self.fig.add_axes([panel_left, y_pos, panel_width, slider_height])
        self.slider_planet_mass = Slider(ax_planet_mass, "Planet Mass", 0.1, 20, 
                                         valinit=self.planet_mass, color="cyan", valstep=0.1, valfmt="%.1f")
        self.slider_planet_mass.on_changed(self.update_planet_mass)
        
        y_pos -= spacing
        ax_planet_dist = self.fig.add_axes([panel_left, y_pos, panel_width, slider_height])
        self.slider_planet_dist = Slider(ax_planet_dist, "Planet Dist", 0.5, 5, 
                                         valinit=self.planet_distance, color="cyan", valstep=0.1, valfmt="%.1f")
        self.slider_planet_dist.on_changed(self.update_planet_dist)
        
        y_pos -= spacing
        ax_planet_vel = self.fig.add_axes([panel_left, y_pos, panel_width, slider_height])
        self.slider_planet_vel = Slider(ax_planet_vel, "Planet Vel", 5, 40, 
                                        valinit=self.planet_velocity, color="cyan", valstep=0.5, valfmt="%.1f")
        self.slider_planet_vel.on_changed(self.update_planet_vel)
        
        y_pos -= spacing * 1.5
        ax_moon_mass = self.fig.add_axes([panel_left, y_pos, panel_width, slider_height])
        self.slider_moon_mass = Slider(ax_moon_mass, "Moon Mass", 0.01, 5, 
                                       valinit=self.moon_mass, color="gray", valstep=0.01, valfmt="%.2f")
        self.slider_moon_mass.on_changed(self.update_moon_mass)
        
        y_pos -= spacing
        ax_moon_dist = self.fig.add_axes([panel_left, y_pos, panel_width, slider_height])
        self.slider_moon_dist = Slider(ax_moon_dist, "Moon Dist", 0.05, 1, 
                                       valinit=self.moon_distance, color="gray", valstep=0.01, valfmt="%.2f")
        self.slider_moon_dist.on_changed(self.update_moon_dist)
        
        y_pos -= spacing
        ax_moon_vel = self.fig.add_axes([panel_left, y_pos, panel_width, slider_height])
        self.slider_moon_vel = Slider(ax_moon_vel, "Moon Vel", 1, 15, 
                                      valinit=self.moon_velocity, color="gray", valstep=0.1, valfmt="%.1f")
        self.slider_moon_vel.on_changed(self.update_moon_vel)
        
        y_pos -= spacing * 1.5
        ax_dt = self.fig.add_axes([panel_left, y_pos, panel_width, slider_height])
        self.slider_dt = Slider(ax_dt, "Time Step", 0.0001, 0.01, 
                               valinit=0.001, color="white", valstep=0.0001, valfmt="%.4f")
        self.slider_dt.on_changed(self.update_dt)
        
        # Buttons
        y_pos -= spacing * 2
        ax_run = self.fig.add_axes([panel_left, y_pos, panel_width/2 - 0.01, 0.04])
        self.btn_run = Button(ax_run, "RUN", color="green", hovercolor="lightgreen")
        self.btn_run.on_clicked(self.toggle_run)
        
        ax_reset = self.fig.add_axes([panel_left + panel_width/2 + 0.01, y_pos, 
                                      panel_width/2 - 0.01, 0.04])
        self.btn_reset = Button(ax_reset, "RESET", color="red", hovercolor="lightcoral")
        self.btn_reset.on_clicked(self.reset_clicked)
        
        # Info Text
        y_pos -= spacing * 2
        info_text = (
            "Key:\n"
            "• Red line = Perfect circular orbit\n"
            "• Blue line = Hill sphere limit (max dist)\n"
            "  (keep moon distance left of the blue line)\n"
            "• Good for blue line to not be there\n\n"
        )
        self.fig.text(panel_left, y_pos - 0.15, info_text, 
                     fontsize=9, color="white", 
                     verticalalignment="top",
                     bbox=dict(boxstyle="round", facecolor="#2a2a2a", alpha=0.8))

        self.update_hints()

    def update_hints(self): # Update Red (Orbit) and Blue (Hill Limit) lines
        
        # planet hints
        try:
            total_mass_sun_sys = self.sun_mass + self.planet_mass
            target_v = np.sqrt(G * total_mass_sun_sys / self.planet_distance)
        except: target_v = 0
            
        try:
            target_r = (G * total_mass_sun_sys) / (self.planet_velocity**2) if self.planet_velocity > 0 else 0
        except: target_r = 0

        # moon hints
        # 1. Circular orbit velocity
        try:
            # For heavy moons, we orbit the Barycenter, so mass is M1+M2
            total_mass_planet_sys = self.planet_mass + self.moon_mass
            target_moon_v = np.sqrt(G * total_mass_planet_sys / self.moon_distance)
        except: target_moon_v = 0
        
        # 2. Hill sphere limit (max dist before stolen by sun)
        # Formula: r_hill = r_planet * cbrt(M_planet / 3*M_sun)
        try:
            # Use planet + moon mass for stronger gravity well
            hill_radius = self.planet_distance * np.cbrt((self.planet_mass + self.moon_mass) / (3 * self.sun_mass))
        except: hill_radius = 0

        # Draw lines helper
        def draw_hint(name, slider, value, color="red"):
            if name in self.hint_lines and self.hint_lines[name]:
                self.hint_lines[name].set_xdata([value, value])
            else:
                line = slider.ax.axvline(value, color=color, linewidth=2, linestyle='--', alpha=0.9)
                self.hint_lines[name] = line

        draw_hint("vel", self.slider_planet_vel, target_v, "red")
        draw_hint("dist", self.slider_planet_dist, target_r, "red")
        
        draw_hint("moon_vel", self.slider_moon_vel, target_moon_v, "red")
        
        # Draw blue hill limit on moon dist
        draw_hint("hill_limit", self.slider_moon_dist, hill_radius, "cyan")
        
        # Red orbit hint on moon dist (standard orbit requirement)
        try:
            target_moon_r = (G * (self.planet_mass + self.moon_mass)) / (self.moon_velocity**2) if self.moon_velocity > 0 else 0
        except: target_moon_r = 0
        draw_hint("moon_dist_target", self.slider_moon_dist, target_moon_r, "red")

        self.fig.canvas.draw_idle()

    # Callbacks
    def update_sun_mass(self, val): self.sun_mass = val; self.update_hints()
    def update_planet_mass(self, val): self.planet_mass = val; self.update_hints()
    def update_planet_dist(self, val): self.planet_distance = val; self.update_hints()
    def update_planet_vel(self, val): self.planet_velocity = val; self.update_hints()
    def update_moon_mass(self, val): self.moon_mass = val; self.update_hints()
    def update_moon_dist(self, val): self.moon_distance = val; self.update_hints()
    def update_moon_vel(self, val): self.moon_velocity = val; self.update_hints()
    def update_dt(self, val): self.dt = val
    
    def toggle_run(self, event):
        self.running = not self.running
        self.btn_run.label.set_text("PAUSE" if self.running else "RUN")
        self.btn_run.color = "orange" if self.running else "green"
    
    def reset_clicked(self, event):
        self.running = False
        self.btn_run.label.set_text("RUN")
        self.btn_run.color = "green"
        self.reset_simulation()
    
    def reset_simulation(self):
        sun = Body("Sun", self.sun_mass, [0, 0, 0], [0, 0, 0], "yellow", 200)
        planet = Body("Planet", self.planet_mass,
                     [self.planet_distance, 0, 0],
                     [0, self.planet_velocity, 0],
                     "cyan", 100)
        moon = Body("Moon", self.moon_mass,
                   [self.planet_distance + self.moon_distance, 0, 0],
                   [0, self.planet_velocity + self.moon_velocity, 0],
                   "lightgray", 50)
        
        self.bodies = [sun, planet, moon]
        
        self.ax.clear()
        self.ax.set_facecolor("black")
        self.ax.set_xlabel("X"); self.ax.set_ylabel("Y"); self.ax.set_zlabel("Z")
        self.ax.tick_params(colors="white")
        
        max_dist = self.planet_distance * 1.5
        self.ax.set_xlim([-max_dist, max_dist])
        self.ax.set_ylim([-max_dist, max_dist])
        self.ax.set_zlim([-max_dist, max_dist])
    
    def animate(self, frame):
        if self.running:
            for _ in range(5): update_system(self.bodies, self.dt)
        
        self.ax.clear()
        self.ax.set_facecolor("black")
        self.ax.set_title(f"3-Body Simulator", color="white")
        
        for body in self.bodies:
            if len(body.trail_x) > 2:
                self.ax.plot(body.trail_x, body.trail_y, body.trail_z,
                           color=body.color, linewidth=1, alpha=0.6)
            self.ax.scatter(body.position[0], body.position[1], body.position[2],
                          color=body.color, s=body.size)
        
        if self.running and self.bodies[1].trail_x:
            # follow planet roughly
            x, y = self.bodies[1].position[0], self.bodies[1].position[1]
            span = self.planet_distance * 0.8
            self.ax.set_xlim([x - span, x + span])
            self.ax.set_ylim([y - span, y + span])
            self.ax.set_zlim([-span, span])
        
        return []
    
    def show(self): plt.show()

if __name__ == "__main__":
    OrbitSimulator().show()