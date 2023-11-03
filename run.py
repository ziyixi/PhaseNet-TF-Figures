from phasenettf.scripts import *
import sys

scripts_mapper = {
    "station_distribution": station_distribution_main,
    "map_with_waveform": map_with_waveform_main,
    "associated_events_horizontal": associated_events_horizontal_main,
    "associated_events_vertical": associated_events_vertical_main,
    "relocated_events_horizontal": relocated_events_horizontal_main,
    "relocated_events_vertical": relocated_events_vertical_main,
    "mantle_wedge_events": mantle_wedge_events_main,
    "manuscript_event_station_distribution": manuscript_event_station_distribution_main,
    "manuscript_workflow_component": manuscript_workflow_component_main,
    "manuscript_phasenet_tf_example": manuscript_phasenet_tf_example_main,
    "manuscript_horizontal_benchmark_1hourwindow": manuscript_horizontal_benchmark_1hourwindow_main,
    "manuscript_horizontal_continious": manuscript_horizontal_continious_main,
    "manuscript_vertical_cross_section_1hour_reference": manuscript_vertical_cross_section_1hour_reference_main,
    "manuscript_vertical_cross_section_1hour_associated": manuscript_vertical_cross_section_1hour_associated_main,
    "manuscript_vertical_cross_section_1hour_relocated": manuscript_vertical_cross_section_1hour_relocated_main,
    "manuscript_vertical_cross_section_1hour_bootstrapped": manuscript_vertical_cross_section_1hour_bootstrapped_main,
    "manuscript_vertical_cross_section_continious_associated": manuscript_vertical_cross_section_continious_associated_main,
    "manuscript_vertical_cross_section_continious_relocated": manuscript_vertical_cross_section_continious_relocated_main,
    "manuscript_vertical_cross_section_continious_bootstrapped": manuscript_vertical_cross_section_continious_bootstrapped_main,
    "manuscript_vertical_cross_section_continious_semi": manuscript_vertical_cross_section_continious_semi_main,
    "manuscript_histogram_1hour_catalog": manuscript_histogram_1hour_catalog_main,
    "manuscript_histogram_1hour_phase": manuscript_histogram_1hour_phase_main,
    "manuscript_histogram_continious_phase": manuscript_histogram_continious_phase_main,
    "manuscript_histogram_continious_catalog": manuscript_histogram_continious_catalog_main,
    "ptf_example_with_ps": manuscript_ptf_example_with_ps_main,
    "horizontal_top100km": horizontal_top100km_main,
    "vertical_cross_section_zoom": vertical_cross_section_zoom_main,
    "manuscript_gamma1d_example": manuscript_gamma1d_example_main,
    "manuscript_continuous_waveform": manuscript_continuous_waveform_main,
    "manuscript_compare_models_residuals": manuscript_compare_models_residuals_main,
    "manuscript_compare_gamma_gamma1d_residuals": manuscript_compare_gamma_gamma1d_residuals_main,
    "manuscript_vertical_crosssections_0to700km": manuscript_vertical_crosssections_0to700km_main,
    "manuscript_vertical_crosssections_0to300km": manuscript_vertical_crosssections_0to300km_main,
    "manuscript_vertical_crosssections_300to700km": manuscript_vertical_crosssections_300to700km_main,
    "manuscript_horizontal_continious_300to400km": manuscript_horizontal_continious_300to400km_main,
    "manuscirpt_final_catalog_plot": manuscirpt_final_catalog_plot_main,
    "manuscirpt_final_catalog_ehb_plot": manuscirpt_final_catalog_ehb_plot_main,
    "manuscript_vertical_continuous": manuscript_vertical_continuous_main,
    "manuscript_compare_final_ehb_catalogs": manuscript_compare_final_ehb_catalogs_main,
    "manualscript_final_catalog_three_columns_version": manualscript_final_catalog_three_columns_version_main,
}


def main():
    if len(sys.argv) == 2:
        if sys.argv[1] == "all" or (sys.argv[1] in scripts_mapper):
            if sys.argv[1] != "all":
                scripts_mapper[sys.argv[1]]()
            else:
                for key in scripts_mapper:
                    print(f"Plot {key} now...")
                    scripts_mapper[key]()
        else:
            raise Exception(f"scripts {sys.argv[1]} is not supported!")
    else:
        raise Exception("correct format: python run.py [script name]")


if __name__ == "__main__":
    main()
