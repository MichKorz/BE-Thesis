import argparse
import os
import matplotlib.pyplot as plt
import LightwaveExplorer as lwe
import numpy as np
from scipy.ndimage import gaussian_filter1d
from colors import colors

def get_palette_color(color_string):
    if color_string.endswith(':alt'):
        base_name = color_string.replace(':alt', '')
        try:
            return getattr(colors, base_name).alternative
        except AttributeError:
            return color_string
    else:
        try:
            return getattr(colors, color_string).base
        except AttributeError:
            return color_string

def main():
    parser = argparse.ArgumentParser(description="Plot Lightwave Explorer simulations.")

    parser.add_argument('-i', '--inputs', nargs='+', required=True,
                        help="Simulation files and optional import limits (e.g., -i sim1.zip 250 2000 sim2.zip)")
    parser.add_argument('-c', '--colors', action='append', nargs='+',
                        help="Syntax: -c [source_index] [base_color] [start1] [end1] [color1] ... (e.g., -c 0 blue:alt 250 500 red:alt)")
    parser.add_argument('-t', '--thickness', action='append', nargs=2,
                        help="Syntax: -t [source_index] [thickness] (e.g., -t 0 2.5). Default thickness is 1.0.")
    parser.add_argument('-u', '--unit', choices=['freq', 'wlen'], default='freq', help="Base unit: freq (THz) or wlen (nm)")
    parser.add_argument('-s', '--smoothness', type=float, default=0.0, help="Sigma for Gaussian smoothing")
    parser.add_argument('-o', '--output', type=str, default=None, help="Output image filename")
    parser.add_argument('-x', '--xlim', nargs=2, type=float, help="Set x-axis limits (e.g., -x 400 1200)")
    parser.add_argument('-n', '--normalize', nargs='*', type=float, help="Normalize globally, or boost specific regions to 1 while scaling the rest independently (e.g., -n 250 500)")
    parser.add_argument('-ne', '--normalize_exclusive', nargs='+', type=float, help="Normalize specific regions to 1 and hide everything outside them (e.g., -ne 250 500)")
    parser.add_argument('-m', '--xticks', nargs='+', type=float, help="Specify custom x-axis tick marks (e.g., -m 400 800 1000 1200)")
    parser.add_argument('-al', '--axis_linewidth', type=float, default=0.8,
                        help="Thickness of axis lines and tick marks. Default is 0.8.")
    parser.add_argument('-ls', '--label_size', nargs='+', default=[10.0, 'normal'],
                        help="Font size and optional text thickness (e.g., -ls 12 bold). Defaults: 10.0 normal.")
    parser.add_argument('-ts', '--tick_size', nargs='+', default=[10.0, 'normal'],
                        help="Font size and optional text thickness (e.g., -ts 10 bold). Defaults: 10.0 normal.")
    parser.add_argument('-d', '--dump_csv', action='store_true',
                        help="Dump the processed and bounded x and y data to a CSV file (semicolon separated, no header) before plotting.")
    parser.add_argument('-in', '--integrate', action='append', nargs=3,
                        help="Syntax: -in [source_index] [start] [end] (e.g., -in 0 400 600). Integrates the area under the curve for the specified source over the given range, shades it, and adds the value to the legend.")

    args = parser.parse_args()

    try:
        ls_val = float(args.label_size[0])
    except ValueError:
        ls_val = 10.0
    ls_weight = args.label_size[1] if len(args.label_size) > 1 else 'normal'

    try:
        ts_val = float(args.tick_size[0])
    except ValueError:
        ts_val = 10.0
    ts_weight = args.tick_size[1] if len(args.tick_size) > 1 else 'normal'

    plt.rcParams['axes.linewidth'] = args.axis_linewidth
    plt.rcParams['xtick.major.width'] = args.axis_linewidth
    plt.rcParams['ytick.major.width'] = args.axis_linewidth
    plt.rcParams['xtick.minor.width'] = args.axis_linewidth
    plt.rcParams['ytick.minor.width'] = args.axis_linewidth
    plt.rcParams['axes.labelsize'] = ls_val
    plt.rcParams['axes.labelweight'] = ls_weight
    plt.rcParams['axes.titlesize'] = ls_val
    plt.rcParams['axes.titleweight'] = ls_weight
    plt.rcParams['xtick.labelsize'] = ts_val
    plt.rcParams['ytick.labelsize'] = ts_val

    plot_xmin = None
    plot_xmax = None
    if args.xlim:
        plot_xmin, plot_xmax = min(args.xlim), max(args.xlim)
    elif args.unit == 'freq':
        plot_xmin, plot_xmax = 0, 850

    parsed_inputs = []
    i = 0
    while i < len(args.inputs):
        filename = args.inputs[i]
        min_val = None
        max_val = None

        if i + 2 < len(args.inputs):
            try:
                min_val = float(args.inputs[i+1])
                max_val = float(args.inputs[i+2])
                i += 3
            except ValueError:
                i += 1
        else:
            i += 1

        parsed_inputs.append((filename, min_val, max_val))

    color_configs = {}
    if args.colors:
        for c_arg in args.colors:
            try:
                src_idx = int(c_arg[0])
                base_col = c_arg[1]
                ranges = []

                for i_c in range(2, len(c_arg), 3):
                    if i_c + 2 < len(c_arg):
                        ranges.append((float(c_arg[i_c]), float(c_arg[i_c+1]), c_arg[i_c+2]))
                    else:
                        pass

                color_configs[src_idx] = {'base': base_col, 'ranges': ranges}
            except (ValueError, IndexError):
                pass

    thickness_configs = {}
    if args.thickness:
        for t_arg in args.thickness:
            try:
                src_idx = int(t_arg[0])
                lw = float(t_arg[1])
                thickness_configs[src_idx] = lw
            except ValueError:
                pass

    integrate_configs = {}
    if args.integrate:
        for in_arg in args.integrate:
            try:
                src_idx = int(in_arg[0])
                start_val = float(in_arg[1])
                end_val = float(in_arg[2])
                if src_idx not in integrate_configs:
                    integrate_configs[src_idx] = []
                integrate_configs[src_idx].append((min(start_val, end_val), max(start_val, end_val)))
            except ValueError:
                pass

    plt.figure(figsize=(12, 4))

    y_axis_label = "Intensity"

    for idx, (sim_file, min_val, max_val) in enumerate(parsed_inputs):
        if not os.path.exists(sim_file):
            continue

        results = lwe.load(sim_file, False)
        y_data = results.spectrumTotal

        if args.smoothness > 0:
            y_data = gaussian_filter1d(y_data, sigma=args.smoothness)

        config = color_configs.get(idx, {'base': 'blue', 'ranges': []})
        base_color = get_palette_color(config['base'])
        line_width = thickness_configs.get(idx, 1.0)

        if args.unit == 'wlen':
            c = 299792458.0
            freq = results.frequencyVectorSpectrum
            valid_idx = freq > 0
            freq_valid = freq[valid_idx]

            x_data = (c / freq_valid) * 1e9
            y_data = y_data[valid_idx] * (freq_valid**2 / c)

            xlabel = "Wavelength (nm)"
            ylabel_base = "Intensity (J/m)"
        else:
            x_data = results.frequencyVectorSpectrum / 1e12
            xlabel = "Frequency (THz)"
            ylabel_base = "Intensity (J/THz)"

        if min_val is not None and max_val is not None:
            start_bound = min(min_val, max_val)
            end_bound = max(min_val, max_val)
            bound_mask = (x_data >= start_bound) & (x_data <= end_bound)
            x_data = x_data[bound_mask]
            y_data = y_data[bound_mask]

            if len(x_data) == 0:
                continue

        if args.normalize_exclusive is not None:
            y_axis_label = "Normalized Intensity (a.u.)"

            new_y_data = np.full_like(y_data, np.nan, dtype=float)

            for i_ne in range(0, len(args.normalize_exclusive) - 1, 2):
                start = min(args.normalize_exclusive[i_ne], args.normalize_exclusive[i_ne+1])
                end = max(args.normalize_exclusive[i_ne], args.normalize_exclusive[i_ne+1])

                mask = (x_data >= start) & (x_data <= end)
                if np.any(mask):
                    local_max = np.nanmax(y_data[mask])
                    if local_max > 0:
                        new_y_data[mask] = y_data[mask] / local_max
                    else:
                        new_y_data[mask] = y_data[mask]

            y_data = new_y_data

        elif args.normalize is not None:
            y_axis_label = "Normalized Intensity (a.u.)"

            if len(args.normalize) == 0:
                global_max = np.nanmax(y_data)
                if global_max > 0:
                    y_data = y_data / global_max
            else:
                inside_mask = np.zeros_like(x_data, dtype=bool)

                for i_n in range(0, len(args.normalize) - 1, 2):
                    start = min(args.normalize[i_n], args.normalize[i_n+1])
                    end = max(args.normalize[i_n], args.normalize[i_n+1])

                    mask = (x_data >= start) & (x_data <= end)
                    if np.any(mask):
                        inside_mask |= mask
                        local_max = np.nanmax(y_data[mask])
                        if local_max > 0:
                            y_data[mask] = y_data[mask] / local_max

                outside_mask = ~inside_mask
                if np.any(outside_mask):
                    outside_max = np.nanmax(y_data[outside_mask])
                    if outside_max > 0:
                        y_data[outside_mask] = y_data[outside_mask] / outside_max
        else:
            y_axis_label = ylabel_base

        if args.dump_csv:
            dump_x = x_data
            dump_y = y_data

            if plot_xmin is not None and plot_xmax is not None:
                x_mask = (dump_x >= plot_xmin) & (dump_x <= plot_xmax)
                dump_x = dump_x[x_mask]
                dump_y = dump_y[x_mask]

            valid_mask = ~np.isnan(dump_y)
            dump_x = dump_x[valid_mask]
            dump_y = dump_y[valid_mask]

            base_filename = os.path.splitext(os.path.basename(sim_file))[0]
            dump_filename = f"{base_filename}_plot_data.csv"
            np.savetxt(dump_filename, np.column_stack((dump_x, dump_y)), delimiter=';', fmt='%g')

        base_label = os.path.basename(sim_file)
        if idx in integrate_configs:
            for start_val, end_val in integrate_configs[idx]:
                mask = (x_data >= start_val) & (x_data <= end_val)
                valid_mask = ~np.isnan(y_data)
                combined_mask = mask & valid_mask

                if np.any(combined_mask):
                    x_shade = x_data[combined_mask]
                    y_shade = y_data[combined_mask]

                    sort_idx = np.argsort(x_shade)
                    x_shade = x_shade[sort_idx]
                    y_shade = y_shade[sort_idx]

                    area = np.trapezoid(y_shade, x_shade)
                    plt.fill_between(x_shade, 0, y_shade, color=base_color, alpha=0.3)

                    base_label += f" | Area: {area:.3g}"

        plt.plot(x_data, y_data, color=base_color, label=base_label, linewidth=line_width, linestyle='-')

        for r_start, r_end, r_col in config['ranges']:
            overlay_color = get_palette_color(r_col)

            mask = (x_data >= r_start) & (x_data <= r_end)
            if np.any(mask):
                plt.plot(x_data[mask], y_data[mask], color=overlay_color, linewidth=line_width, linestyle='-')

    plt.xlabel(xlabel)
    plt.ylabel(y_axis_label)

    if plot_xmin is not None and plot_xmax is not None:
        plt.xlim(plot_xmin, plot_xmax)

    if args.xticks:
        plt.xticks(args.xticks)

    plt.title("Simulation Results")
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.3)

    ax = plt.gca()
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight(ts_weight)

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

if __name__ == "__main__":
    main()
