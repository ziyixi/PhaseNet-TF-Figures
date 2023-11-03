from .station_distribution import main as station_distribution_main
from .map_with_waveform import main as map_with_waveform_main
from .associated_events_horizontal import main as associated_events_horizontal_main
from .associated_events_vertical import main as associated_events_vertical_main
from .relocated_events_horizontal import main as relocated_events_horizontal_main
from .relocated_events_vertical import main as relocated_events_vertical_main
from .mantle_wedge_events import main as mantle_wedge_events_main
from .manuscript_event_station_distribution import (
    main as manuscript_event_station_distribution_main,
)
from .manuscript_workflow_component import main as manuscript_workflow_component_main
from .manuscript_phasenet_tf_example import main as manuscript_phasenet_tf_example_main
from .manuscript_horizontal_benchmark_1hourwindow import (
    main as manuscript_horizontal_benchmark_1hourwindow_main,
)
from .manuscript_horizontal_continious import (
    main as manuscript_horizontal_continious_main,
)
from .manuscript_vertical_cross_section_1hour_reference import (
    main as manuscript_vertical_cross_section_1hour_reference_main,
)
from .manuscript_vertical_cross_section_1hour_associated import (
    main as manuscript_vertical_cross_section_1hour_associated_main,
)
from .manuscript_vertical_cross_section_1hour_relocated import (
    main as manuscript_vertical_cross_section_1hour_relocated_main,
)
from .manuscript_vertical_cross_section_1hour_bootstrapped import (
    main as manuscript_vertical_cross_section_1hour_bootstrapped_main,
)
from .manuscript_vertical_cross_section_continious_associated import (
    main as manuscript_vertical_cross_section_continious_associated_main,
)
from .manuscript_vertical_cross_section_continious_relocated import (
    main as manuscript_vertical_cross_section_continious_relocated_main,
)
from .manuscript_vertical_cross_section_continious_bootstrapped import (
    main as manuscript_vertical_cross_section_continious_bootstrapped_main,
)
from .manuscript_vertical_cross_section_continious_semi import (
    main as manuscript_vertical_cross_section_continious_semi_main,
)
from .manuscript_histogram_1hour_catalog import (
    main as manuscript_histogram_1hour_catalog_main,
)
from .manuscript_histogram_1hour_phase import (
    main as manuscript_histogram_1hour_phase_main,
)
from .manuscript_histogram_continious_phase import (
    main as manuscript_histogram_continious_phase_main,
)
from .manuscript_histogram_continious_catalog import (
    main as manuscript_histogram_continious_catalog_main,
)
from .ptf_example_with_ps import main as manuscript_ptf_example_with_ps_main
from .horizontal_top100km import main as horizontal_top100km_main
from .vertical_cross_section_zoom import main as vertical_cross_section_zoom_main
from .manuscript_gamma1d_example import main as manuscript_gamma1d_example_main
from .manuscript_continuous_waveform import main as manuscript_continuous_waveform_main
from .manuscript_compare_models_residuals import (
    main as manuscript_compare_models_residuals_main,
)
from .manuscript_compare_gamma_gamma1d_residuals import (
    main as manuscript_compare_gamma_gamma1d_residuals_main,
)
from .manuscript_vertical_crosssections_0to700km import (
    main as manuscript_vertical_crosssections_0to700km_main,
)
from .manuscript_vertical_crosssections_0to300km import (
    main as manuscript_vertical_crosssections_0to300km_main,
)
from .manuscript_vertical_crosssections_300to700km import (
    main as manuscript_vertical_crosssections_300to700km_main,
)
from .manuscript_horizontal_continious_300to400km import (
    main as manuscript_horizontal_continious_300to400km_main,
)
from .manuscirpt_final_catalog_plot import main as manuscirpt_final_catalog_plot_main
from .manuscirpt_final_catalog_ehb_plot import (
    main as manuscirpt_final_catalog_ehb_plot_main,
)
from .manuscript_vertical_continuous import main as manuscript_vertical_continuous_main
from .manuscript_compare_final_ehb_catalogs import (
    main as manuscript_compare_final_ehb_catalogs_main,
)
from .manualscript_final_catalog_three_columns_version import (
    main as manualscript_final_catalog_three_columns_version_main,
)

__all__ = [
    "station_distribution_main",
    "map_with_waveform_main",
    "associated_events_horizontal_main",
    "associated_events_vertical_main",
    "relocated_events_horizontal_main",
    "relocated_events_vertical_main",
    "mantle_wedge_events_main",
    "manuscript_event_station_distribution_main",
    "manuscript_workflow_component_main",
    "manuscript_phasenet_tf_example_main",
    "manuscript_horizontal_benchmark_1hourwindow_main",
    "manuscript_horizontal_continious_main",
    "manuscript_vertical_cross_section_1hour_reference_main",
    "manuscript_vertical_cross_section_1hour_associated_main",
    "manuscript_vertical_cross_section_1hour_relocated_main",
    "manuscript_vertical_cross_section_1hour_bootstrapped_main",
    "manuscript_vertical_cross_section_continious_associated_main",
    "manuscript_vertical_cross_section_continious_relocated_main",
    "manuscript_vertical_cross_section_continious_bootstrapped_main",
    "manuscript_vertical_cross_section_continious_semi_main",
    "manuscript_histogram_1hour_catalog_main",
    "manuscript_histogram_1hour_phase_main",
    "manuscript_histogram_continious_phase_main",
    "manuscript_histogram_continious_catalog_main",
    "manuscript_ptf_example_with_ps_main",
    "horizontal_top100km_main",
    "vertical_cross_section_zoom_main",
    "manuscript_gamma1d_example_main",
    "manuscript_continuous_waveform_main",
    "manuscript_compare_models_residuals_main",
    "manuscript_compare_gamma_gamma1d_residuals_main",
    "manuscript_vertical_crosssections_0to700km_main",
    "manuscript_vertical_crosssections_0to300km_main",
    "manuscript_vertical_crosssections_300to700km_main",
    "manuscript_horizontal_continious_300to400km_main",
    "manuscirpt_final_catalog_plot_main",
    "manuscirpt_final_catalog_ehb_plot_main",
    "manuscript_vertical_continuous_main",
    "manuscript_compare_final_ehb_catalogs_main",
    "manualscript_final_catalog_three_columns_version_main",
]
