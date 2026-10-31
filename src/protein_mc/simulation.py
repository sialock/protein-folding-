"""Single-chain Metropolis Monte Carlo with independent seeded RNG state."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .energy import GoModel
from .geometry import crankshaft_move, kabsch_rmsd


@dataclass(frozen=True)
class RunConfig:
    steps: int = 1500
    scramble_steps: int = 300
    temperature: float = 1.0  # dimensionless, NOT Kelvin
    max_angle_degrees: float = 40.0
    seed: int = 42
    report_every: int = 25

    def __post_init__(self) -> None:
        if self.steps < 1 or self.scramble_steps < 0 or self.report_every < 1:
            raise ValueError("steps >= 1, scramble_steps >= 0, report_every >= 1 required")
        if not np.isfinite(self.temperature) or self.temperature <= 0:
            raise ValueError("temperature must be finite and > 0 (reduced units)")
        if not np.isfinite(self.max_angle_degrees) or not 0 < self.max_angle_degrees <= 180:
            raise ValueError("max_angle_degrees must be between 0 and 180")
        if self.seed < 0:
            raise ValueError("seed must be nonnegative")


@dataclass
class RunResult:
    initial_coords: np.ndarray
    final_coords: np.ndarray
    best_coords: np.ndarray
    trace: list[dict]
    accepted: int
    attempted: int
    best_energy: float

    @property
    def acceptance_rate(self) -> float:
        return self.accepted / self.attempted


def simulate(model: GoModel, config: RunConfig) -> RunResult:
    rng = np.random.default_rng(config.seed)
    coords = model.reference.copy()
    for _ in range(config.scramble_steps):
        coords = crankshaft_move(coords, rng, max_angle_radians=np.deg2rad(110))
    initial = coords.copy()
    energy = model.energy(coords)
    best = coords.copy()
    best_energy = energy
    accepted = 0
    trace: list[dict] = []

    def record(step: int, accepted_move: bool) -> None:
        trace.append(
            {
                "step": step,
                "energy_reduced": energy,
                "rmsd_angstrom": kabsch_rmsd(coords, model.reference),
                "native_contact_fraction": model.native_contact_fraction(coords),
                "accepted_last_move": accepted_move,
                "cumulative_acceptance_rate": accepted / step if step else 0.0,
            }
        )

    record(0, False)
    for step in range(1, config.steps + 1):
        proposal = crankshaft_move(
            coords, rng, max_angle_radians=np.deg2rad(config.max_angle_degrees)
        )
        next_energy = model.energy(proposal)
        delta = next_energy - energy
        did_accept = bool(delta <= 0 or rng.random() < np.exp(-delta / config.temperature))
        if did_accept:
            coords = proposal
            energy = next_energy
            accepted += 1
            if energy < best_energy:
                best = coords.copy()
                best_energy = energy
        if step % config.report_every == 0 or step == config.steps:
            record(step, did_accept)

    return RunResult(initial, coords.copy(), best, trace, accepted, config.steps, best_energy)