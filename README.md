# Parabolic Vessel — Hard Sphere Gas Simulation

A real-time 2D molecular dynamics simulation of hard-sphere particles in a
configurable vessel, written in C with SDL2. Demonstrates thermalisation,
energy conservation, Maxwell-Boltzmann statistics, and acoustic compression
waves via an oscillating piston ceiling.

---

## Physics

Particles are point masses under constant downward gravity inside a vessel
whose floor is a parabolic curve:

```
y_floor(x) = PARA_BASE - para_a * (x - sim_w/2)²
```

Setting `para_a = 0` gives a flat rectangular box. Vertical side walls and a
movable ceiling piston complete the enclosure.

**Collisions are fully elastic** by default (`restitution = 1.0`):

- **Particle–parabola / walls**: continuous collision detection (CCD) via
  binary search finds the exact contact time; velocity is reflected about the
  surface normal.
- **Particle–particle**: fully event-driven MD — all pairwise contact times
  are solved analytically (the quadratic formula; gravity cancels in the
  relative frame so it is exact). The earliest event wins, all particles
  advance to that moment, and the pair exchanges normal-direction velocity
  components.
- **Ceiling piston**: binary search finds exact contact time with the moving
  surface; the bounce reflects velocity relative to the piston, transferring
  the correct momentum.

This approach gives **< 0.02% energy drift over 10 000 frames** with no
piston, compared to ~35% drift from naive post-step overlap detection.

### Thermalisation

A single particle bouncing in the parabola conserves its own energy
indefinitely — the bowl is integrable. Thermalisation only occurs through
**particle–particle collisions**, which redistribute kinetic energy. Starting
monoenergetic (all same speed), the KE distribution broadens and converges to
the 2D Maxwell-Boltzmann exponential:

```
f(KE) ∝ exp(−KE / kT)    where kT = mean KE per particle
```

### Oscillating piston

When `ceil_amp > 0` the ceiling acts as a sinusoidal piston:

```
y_ceil(t) = radius + ceil_amp × (1 + sin(2π × ceil_freq × t))
```

This drives compression waves that propagate visibly through the gas. For
near-zero net energy drift the piston must be **quasi-static**:

```
ceil_amp × 2π × ceil_freq  ≪  mono_speed   (ratio < ~0.1)
```

When the piston moves faster than the gas particles, irreversible **Fermi
acceleration** heats the gas permanently — this is real physics, not a bug.

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

> The source uses `#include <SDL.h>` (flat include). Change to
> `#include <SDL2/SDL.h>` if your setup requires it.

---

## Configuration

Parameters are loaded from a plain-text config file at startup. Most can
also be adjusted live via the **in-app sliders** in the right-hand panel.

### Quick start

```bash
./ballsim --write-config ballsim.cfg   # write a default template
$EDITOR ballsim.cfg                    # edit to taste
./ballsim                              # auto-loads ballsim.cfg
./ballsim --config experiment.cfg      # explicit path
```

### Box dimensions

`sim_w` and `sim_h` **must be set in the config file** — they require a window
resize and cannot be changed via sliders. Edit `ballsim.cfg` and restart.

```ini
sim_w = 900    # simulation area width  (200 – 2560 px)
sim_h = 800    # simulation area height (200 – 1440 px)
```

### Full parameter reference

| Parameter | Default | Slider | Effect |
|-----------|---------|--------|--------|
| `sim_w` | 900 | — | Simulation area width (restart required) |
| `sim_h` | 800 | — | Simulation area height (restart required) |
| `n_particles` | 150 | ✓ | Starting particle count (1–500) |
| `gravity` | 500 px/s² | ✓ | Stronger → faster bouncing, more PE at bottom |
| `radius` | 5 px | ✓ | Larger → more collisions, faster thermalisation |
| `restitution` | 1.0 | ✓ | < 1.0 adds damping; energy trace slopes down |
| `mono_speed` | 220 px/s | ✓ | Speed for M-key monoenergetic reset |
| `para_a` | 0.003 | ✓ | Parabola curvature; 0 = flat rectangular box |
| `ceil_amp` | 0 px | ✓ | Piston amplitude; 0 = fixed ceiling |
| `ceil_freq` | 0.5 Hz | ✓ | Piston oscillation frequency |

Parameters are clamped to safe ranges on load — typos won't crash the sim.
Unknown keys are silently ignored.

---

## In-app sliders

The lower portion of the right-hand panel contains eight sliders covering all
runtime-tunable parameters. Click anywhere on a track to jump to that value;
hold and drag to scrub continuously. The active slider highlights in bright
cyan. Changes that affect particle placement (radius, para_a, n_particles,
mono_speed) trigger an automatic reset.

---

## Controls

| Key | Action |
|-----|--------|
| `M` | Reset — monoenergetic (all same speed, best for thermalisation demo) |
| `R` | Reset — random speeds |
| `+` / `=` | Add a particle (up to 500) |
| `-` | Remove a particle |
| `Q` / `Esc` | Quit |

---

## Stats panel

### Numeric readouts

| Readout | Meaning |
|---------|---------|
| BOX SIZE | Current sim_w × sim_h in pixels |
| PARTICLES | Live particle count |
| TOTAL ENERGY (E/E0) | Energy normalised to initial — 1.0000 = perfectly conserved |
| KINETIC (KE/E) | Kinetic energy as fraction of total |
| POTENTIAL (PE/E) | Potential energy as fraction of total |
| TEMPERATURE (KE/N) | Mean KE per particle — proxy for temperature |
| COLLISIONS | Cumulative particle–particle collision count |

KE/E and PE/E always sum to 1.0. At equilibrium in a parabolic bowl expect
KE/E ≈ 0.55–0.65 depending on how high the particles bounce.

### Energy trace

Ring buffer of the last ~280 samples of E/E0. The white horizontal line marks
1.0 (perfect conservation). Bars colour from green (near 1.0) toward red as
drift grows. The DRIFT label shows the current percentage deviation.

### KE distribution

Histogram of per-particle kinetic energies (cyan bars, exponentially smoothed)
overlaid with the theoretical 2D Maxwell-Boltzmann curve (yellow). Press `M`
to reset to a monoenergetic spike and watch it relax to the exponential shape
over a few seconds as collisions redistribute energy.

---

## Suggested experiments

**Thermalisation** — press `M` for monoenergetic start, watch the KE histogram
spike broaden into the exponential MB curve. The time scale depends on particle
density and radius.

**Flat box** — set `para_a = 0` for a rectangular box with uniform gravity.

**Compression waves** — set `para_a = 0`, `ceil_amp = 5`, `ceil_freq = 0.5`.
Watch density waves propagate downward from the piston on each downstroke.
Keep `ceil_amp × 2π × ceil_freq` well below `mono_speed` to avoid heating.

**Inelastic collapse** — set `restitution = 0.8` and watch kinetic energy
slowly drain into the energy trace. Particles cluster at the bottom.

**Dense gas** — increase `n_particles` and `radius` together; the collision
rate rises sharply and thermalisation becomes almost instantaneous.
