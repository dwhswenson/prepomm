import pandas as pd
import matplotlib.pyplot as plt
from simtk import unit as u

def get_pandas(csv_or_pandas):
    if isinstance(csv_or_pandas, pd.DataFrame):
        return csv_or_pandas
    else:
        return pd.read_csv(csv_or_pandas)

def get_time(df):
    return df['#\"Time (ps)\"']

def plot_pt_equilibration(csv_or_pandas, simulation, pressure_label='volume'):
    box_column = {
        'volume': "Box Volume (nm^3)",
        'density': "Density (g/mL)"
    }[pressure_label]
    temperature_color = 'r'
    box_color = 'b'
    temp_column = "Temperature (K)"
    df = get_pandas(csv_or_pandas)
    time = get_time(df)
    box = df[box_column]
    temp = df[temp_column]
    expected = simulation.integrator.getTemperature().value_in_unit(u.kelvin)
    expected_label = "Expected temperature ({:.1f} K)".format(expected)

    (fig, ax1) = plt.subplots()

    temp_line = ax1.plot(time, temp, color=temperature_color)
    expected_temp_line = ax1.plot(time, [expected] * len(time),
                                  color=temperature_color, ls='--',
                                  label=expected_label)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel(temp_column)
    minimal_ylim = expected - 15, expected + 15
    ylim = ax1.get_ylim()
    new_ylim = [m(minimal_ylim[i], ylim[i]) for (i, m) in enumerate([min, max])]
    ax1.set_ylim(*new_ylim)
    ax1.tick_params('y', colors=temperature_color)

    ax2 = ax1.twinx()
    box_line = ax2.plot(time, box, color=box_color)
    ax2.set_ylabel(box_column)
    ax2.tick_params('y', colors=box_color)

    lines = temp_line + expected_temp_line + box_line
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc=3)

    fig.tight_layout()
    return fig
