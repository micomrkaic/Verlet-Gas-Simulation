# Parabolic Vessel — Hard Sphere Gas Simulation

A real-time 2D molecular dynamics simulation of hard-sphere particles in a
parabolic bowl, written in C with SDL2 for rendering. Demonstrates
thermalisation, energy conservation, and convergence to the Maxwell-Boltzmann
distribution.

---

## Physics

Particles are point masses with radius 5 px moving under constant downward
gravity (500 px/s²). The vessel is a parabolic bowl:

```
y_surface(x) = PARA_BASE - PARA_A * (x - SIM_W/2)²
```

with vertical side walls above the rim and a flat ceiling.

**Collisions are fully elastic** (`RESTITUTION = 1.0`):

- **Particle–wall / parabola**: continuous collision detection (CCD) via binary
  search finds the exact contact time; velocity is reflected about the surface
  normal. No energy is added or removed.
- **Particle–particle**: event-driven MD — at each timestep all pairwise contact
  times are solved analytically (quadratic formula; gravity cancels in the
  relative frame so the result is exact), the earliest event wins, all particles
  are advanced to that moment, then the colliding pair exchanges normal-direction
  velocity components.

This event-driven approach gives **< 0.02% energy drift over 10 000 frames**,
compared to ~35% drift from naive post-step overlap detection.

### Why thermalisation happens

A single particle bouncing in the parabola conserves its own energy
indefinitely — the bowl is an integrable system. Thermalisation only occurs
through **particle–particle collisions**, which redistribute kinetic energy
between particles. Starting from a monoenergetic state (all particles at the
same speed), the KE distribution broadens over a few seconds and converges to
the 2D Maxwell-Boltzmann (exponential) distribution:

```
f(KE) ∝ exp(−KE / kT)    where kT = mean KE per particle
```

---

## Building

### Linux

```bash
gcc -O2 -std=c11 -Wall -Wextra -Wno-unused-parameter \
    -I/usr/include/SDL2 -D_REENTRANT \
    -o ballsim main_collisions.c -lSDL2 -lm
```

### macOS (Homebrew SDL2)

```bash
gcc -O2 -std=c11 \
    -I/opt/homebrew/include/SDL2 -L/opt/homebrew/lib \
    -o ballsim main_collisions.c -lSDL2 -lm
```

Intel Macs: replace `/opt/homebrew` with `/usr/local`.

> The source uses `#include <SDL.h>` (flat include), which is correct for
> Homebrew. Change to `#include <SDL2/SDL.h>` if your setup requires it.

---

## Controls

| Key | Action |
|-----|--------|
| `M` | Reset — monoenergetic (all particles same speed, best for watching thermalisation) |
| `R` | Reset — random speeds |
| `+` / `=` | Add a particle |
| `-` | Remove a particle |
| `Q` / `Esc` | Quit |

---

## Stats Panel

The right-hand panel shows live readouts and charts.

| Readout | Meaning |
|---------|---------|
| PARTICLES | Current particle count |
| TOTAL ENERGY (E/E0) | Current energy normalised to initial energy — should stay at 1.0000 |
| KINETIC (KE/E) | Kinetic energy as a fraction of total — watch this equilibrate |
| POTENTIAL (PE/E) | Potential energy as a fraction of total |
| TEMPERATURE (KE/N) | Mean kinetic energy per particle — proxy for temperature |
| COLLISIONS | Cumulative particle–particle collision count |

**Energy trace** — ring buffer of the last ~280 samples of E/E0. The white
horizontal line marks 1.0 (perfect conservation). Bars are green near 1.0 and
shift toward red as drift grows. A "DRIFT" label shows the current percentage
deviation.

**KE distribution histogram** — 30-bin histogram of per-particle kinetic
energies (cyan bars, exponentially smoothed), overlaid with the theoretical 2D
Maxwell-Boltzmann curve (yellow line). Press `M` to reset to a monoenergetic
spike and watch it relax to the exponential shape.

---

## Tuning

Key `#define` values near the top of `main_collisions.c`:

| Constant | Default | Effect |
|----------|---------|--------|
| `N_PARTICLES` | 150 | Starting particle count |
| `GRAVITY` | 500.0 px/s² | Gravitational acceleration |
| `RADIUS` | 5.0 px | Particle radius |
| `RESTITUTION` | 1.0 | Coefficient of restitution (1 = elastic) |
| `MONO_SPEED` | 220.0 px/s | Speed used for monoenergetic reset |
| `PARA_A` | 0.003 | Parabola curvature (larger = steeper bowl) |
