"""
Astrophysical Benchmarks and Work-Precision Pareto Scaling Package.

Exposes:
- BENCHMARK_REGISTRY, get_benchmark
- create_radial_collision_benchmark, create_burrau_pythagorean_benchmark,
  create_arenstorf_lunar_benchmark, create_figure8_choreography_benchmark,
  create_sitnikov_oscillatory_benchmark, create_quadruple_collision_benchmark,
  create_jovian_solar_system_benchmark, create_pareto_test_benchmark
- BenchmarkResult, BenchmarkRunner, compute_work_precision_pareto
"""

from .registry import (
    BENCHMARK_REGISTRY,
    get_benchmark,
    create_radial_collision_benchmark,
    create_burrau_pythagorean_benchmark,
    create_arenstorf_lunar_benchmark,
    create_figure8_choreography_benchmark,
    create_sitnikov_oscillatory_benchmark,
    create_quadruple_collision_benchmark,
    create_jovian_solar_system_benchmark,
    create_pareto_test_benchmark,
    create_eccentric_binary_benchmark,
    create_jovian_laplace_resonance_benchmark,
    create_comet_lexell_1770_benchmark,
    create_shoemaker_levy_9_1992_benchmark,
    create_apophis_2029_flyby_benchmark,
    create_hierarchical_3d_eccentric_triple_benchmark,
    create_inner_solar_system_mercury_benchmark,
)
from .runner import (
    BenchmarkResult,
    BenchmarkRunner,
    compute_work_precision_pareto,
)

__all__ = [
    "BENCHMARK_REGISTRY",
    "get_benchmark",
    "create_radial_collision_benchmark",
    "create_burrau_pythagorean_benchmark",
    "create_arenstorf_lunar_benchmark",
    "create_figure8_choreography_benchmark",
    "create_sitnikov_oscillatory_benchmark",
    "create_quadruple_collision_benchmark",
    "create_jovian_solar_system_benchmark",
    "create_pareto_test_benchmark",
    "create_eccentric_binary_benchmark",
    "create_jovian_laplace_resonance_benchmark",
    "create_comet_lexell_1770_benchmark",
    "create_shoemaker_levy_9_1992_benchmark",
    "create_apophis_2029_flyby_benchmark",
    "create_hierarchical_3d_eccentric_triple_benchmark",
    "create_inner_solar_system_mercury_benchmark",
    "BenchmarkResult",
    "BenchmarkRunner",
    "compute_work_precision_pareto",
]
