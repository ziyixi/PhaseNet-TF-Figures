from phasenettf.scripts import *
import sys

scripts_mapper = {
    "station_distribution": station_distribution_main,
    "map_with_waveform": map_with_waveform_main,
    "associated_events_horizontal": associated_events_horizontal_main,
    "associated_events_vertical": associated_events_vertical_main,
    "relocated_events_horizontal": relocated_events_horizontal_main,
    "relocated_events_vertical": relocated_events_vertical_main,
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
