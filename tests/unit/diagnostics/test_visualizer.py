"""
Unit Tests for Layer 6: Publication 10-Panel Visualizer in Headless Mode (visualizer.py).

Theoretical Alignments:
- Monograph: Chapter 6 and Appendix C (Visualizing Conormal Desingularization).
- Invariants: Headless generation without X11 server (matplotlib.use('Agg')),
  exactly 10 panels structured in a 5x2 grid, valid PNG export with magic header.
"""

import os
import tempfile
import numpy as np
import pytest
import matplotlib.pyplot as plt

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator
from normal_cone_nbody_float64.diagnostics.visualizer import PublicationVisualizer


def test_publication_visualizer_renders_10_panels_headless(tmp_path):
    """
    Validates rendering a full 10-panel dashboard on a simulated 3-body system:
    - Verifies headless backend is active ('Agg').
    - Verifies exactly 10 subplots (5 rows x 2 columns).
    - Verifies saving to disk produces a valid PNG file (> 10 KB with PNG magic bytes).
    """
    assert plt.get_backend().lower() == 'agg'

    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 5.0, 0.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.6, 0.0],
        [0.0, 0.6, 0.0],
        [0.1, 0.0, 0.0],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)

    system = NBodySystem(state=state0, integrator=StrangIntegrator(G=1.0), G=1.0)
    system.evolve(t_span=0.1, dt=0.02)

    visualizer = PublicationVisualizer(dpi=150)
    output_png = str(tmp_path / "test_10panel_output.png")

    fig = visualizer.render_10panel_dashboard(
        history=system.history,
        title="Unit Test 10-Panel Verification",
        save_path=output_png,
    )

    try:
        assert isinstance(fig, plt.Figure)
        # 10 panels in the grid
        assert len(fig.axes) == 10

        # File validation
        assert os.path.exists(output_png)
        file_size = os.path.getsize(output_png)
        assert file_size > 10000, f"Expected non-trivial PNG size, got {file_size} bytes"

        # Check PNG header magic bytes: 0x89 0x50 0x4E 0x47 0x0D 0x0A 0x1A 0x0A
        with open(output_png, "rb") as f:
            header = f.read(8)
            assert header == b"\x89PNG\r\n\x1a\n", "Generated file is not a valid PNG."
    finally:
        visualizer.close(fig)


def test_publication_visualizer_renders_multi_integrator_dashboard(tmp_path):
    """
    Validates rendering the 6-panel multi-integrator invariants dashboard:
    - Verifies exactly 6 subplots (3 rows x 2 columns).
    - Verifies saving to disk produces a valid PNG file (> 10 KB with PNG magic bytes).
    """
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 5.0, 0.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.6, 0.0],
        [0.0, 0.6, 0.0],
        [0.1, 0.0, 0.0],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)

    system_strang = NBodySystem(state=state0, integrator=StrangIntegrator(G=1.0), G=1.0)
    system_strang.evolve(t_span=0.1, dt=0.02)

    histories = {
        "strang": system_strang.history,
    }

    visualizer = PublicationVisualizer(dpi=150)
    output_png = str(tmp_path / "test_multi_invariants_output.png")

    fig = visualizer.render_multi_integrator_invariants_dashboard(
        histories=histories,
        title="Unit Test Multi-Integrator Invariants Verification",
        save_path=output_png,
    )

    try:
        assert isinstance(fig, plt.Figure)
        assert len(fig.axes) == 6

        assert os.path.exists(output_png)
        file_size = os.path.getsize(output_png)
        assert file_size > 10000, f"Expected non-trivial PNG size, got {file_size} bytes"

        with open(output_png, "rb") as f:
            header = f.read(8)
            assert header == b"\x89PNG\r\n\x1a\n", "Generated file is not a valid PNG."
    finally:
        visualizer.close(fig)


def test_publication_visualizer_renders_refined_dashboard(tmp_path):
    """
    Validates rendering the refined 8-panel master dashboard:
    - Verifies valid PNG creation with correct header magic bytes.
    - Verifies support for NBodySystem instance with desingularization_mode.
    """
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 5.0, 0.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.6, 0.0],
        [0.0, 0.6, 0.0],
        [0.1, 0.0, 0.0],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)

    system = NBodySystem(
        state=state0,
        integrator=StrangIntegrator(G=1.0),
        G=1.0,
        desingularization_mode="auto",
    )
    system.evolve(t_span=0.1, dt=0.02)

    visualizer = PublicationVisualizer(dpi=100)
    output_png = str(tmp_path / "test_refined_dashboard_output.png")

    fig = visualizer.render_refined_dashboard(
        system_or_history=system,
        title="Unit Test Refined Dashboard Verification",
        case_name="Test System",
        save_path=output_png,
    )

    try:
        assert isinstance(fig, plt.Figure)
        assert os.path.exists(output_png)
        file_size = os.path.getsize(output_png)
        assert file_size > 10000, f"Expected non-trivial PNG size, got {file_size} bytes"

        with open(output_png, "rb") as f:
            header = f.read(8)
            assert header == b"\x89PNG\r\n\x1a\n", "Generated file is not a valid PNG."
    finally:
        visualizer.close(fig)


def test_visualizer_purges_obsolete_processed_bm6():
    """
    Validates that visualizer.py no longer references the invalidated 'processed_bm6'
    or 'PBM6' scheme (Theorem 5.14, D-146-05).
    """
    import inspect
    from normal_cone_nbody_float64.diagnostics import visualizer as vis_mod

    src = inspect.getsource(vis_mod)
    assert "processed_bm6" not in src
    assert "PBM6" not in src

