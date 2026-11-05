"""
Layer 5: System Package (`normal_cone_nbody_float64.system`).

Provides high-level system orchestration, dynamical atlas maintenance under gamma = 0.6 hysteresis,
discrete dyadic multi-rate scheduling, and trajectory diagnostic logging:
1. `atlas_manager.py`: Deterministic Jacobi atlas management, admissibility checking, and tree flips.
2. `dyadic_scheduler.py`: Discrete dyadic time lattice, colored window reservations, and dynamic timescale adaptation.
3. `history.py`: Historical snapshot ledger and Lie--Bargmann invariant metric tracking.
4. `nbody_system.py`: Unified high-level simulation facade connecting all layers.
"""

from .atlas_manager import AtlasManager, AtlasTransitionRecord
from .dyadic_scheduler import DyadicScheduler, ColoredInterval
from .history import HistoryTracker, HistorySnapshot
from .nbody_system import NBodySystem

__all__ = [
    "AtlasManager",
    "AtlasTransitionRecord",
    "DyadicScheduler",
    "ColoredInterval",
    "HistoryTracker",
    "HistorySnapshot",
    "NBodySystem",
]
