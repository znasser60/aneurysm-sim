import numpy as np
import matplotlib.pyplot as plt


def plot_pressure_vs_stretch(results): 
    sv_stretch_var = results["sv_stretch_var"]
    sv_pressure_var = results["sv_pressure_var"]
    sv_pressure_var_elastin = results["sv_pressure_var_elastin"]
    sv_pressure_var_collagen = results["sv_pressure_var_collagen"]
    sv_pressure_var_collagen_me = results["sv_pressure_var_collagen_me"]
    sv_pressure_var_collagen_ad = results["sv_pressure_var_collagen_ad"]
    sv_pressure_var_muscle_p = results["sv_pressure_var_muscle_p"]
    sv_pressure_var_muscle_a = results["sv_pressure_var_muscle_a"]

    n_zoom = 120

    plt.figure(figsize=(12, 8))
    plt.plot(sv_stretch_var[:n_zoom], sv_pressure_var[:n_zoom]/1e3, linewidth=2, label='Total')
    plt.plot(sv_stretch_var[32:n_zoom], sv_pressure_var_elastin[32:n_zoom]/1e3, '--', linewidth=2, label='Elastin')
    plt.plot(sv_stretch_var[64:n_zoom], sv_pressure_var_collagen[64:n_zoom]/1e3, '-.', linewidth=2, label='Collagen')
    plt.plot(sv_stretch_var[64:n_zoom], sv_pressure_var_collagen_me[64:n_zoom]/1e3, '--', linewidth=2, label='Collagen Media')
    plt.plot(sv_stretch_var[64:n_zoom], sv_pressure_var_collagen_ad[64:n_zoom]/1e3, '--', linewidth=2, label='Collagen Adventitia')
    plt.plot(sv_stretch_var[44:n_zoom], sv_pressure_var_muscle_p[44:n_zoom]/1e3, '--', linewidth=2, label='Muscle Passive')
    plt.plot(sv_stretch_var[:n_zoom], sv_pressure_var_muscle_a[:n_zoom]/1e3, '--', linewidth=2, label='Muscle Active')
    plt.title('Pressure vs Stretch')
    plt.xlabel('Stretch')
    plt.ylabel('Pressure (kPa)')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.ylim([0, 60])
    plt.show()

def plot_pressure_vs_diameter(results): 
    n_zoom_diam = 125
    
    sv_diam_var = results["sv_diam_var"]
    sv_pressure_var = results["sv_pressure_var"]
    sv_pressure_var_elastin = results["sv_pressure_var_elastin"]
    sv_pressure_var_collagen_me = results["sv_pressure_var_collagen_me"]
    sv_pressure_var_collagen_ad = results["sv_pressure_var_collagen_ad"]
    sv_pressure_var_collagen = results["sv_pressure_var_collagen"]
    sv_pressure_var_muscle_p = results["sv_pressure_var_muscle_p"]
    sv_pressure_var_muscle_a = results["sv_pressure_var_muscle_a"]

    sv_pressure_var_coll = sv_pressure_var_collagen_me + sv_pressure_var_collagen_ad

    plt.figure(figsize=(12, 8))
    plt.plot(sv_diam_var[:n_zoom_diam], sv_pressure_var[:n_zoom_diam]/1e3, 'k-', linewidth=6, label='Total')
    plt.plot(sv_diam_var[32:n_zoom_diam], sv_pressure_var_elastin[32:n_zoom_diam]/1e3, 'k--', linewidth=3, label='E')
    plt.plot(sv_diam_var[66:n_zoom_diam], sv_pressure_var_coll[66:n_zoom_diam]/1e3, 'k-.', linewidth=3, label='C')
    plt.plot(sv_diam_var[44:n_zoom_diam], sv_pressure_var_muscle_p[44:n_zoom_diam]/1e3, 'k+', linewidth=2, markersize=8, label='VSMCp')
    plt.plot(sv_diam_var[:n_zoom_diam], sv_pressure_var_muscle_a[:n_zoom_diam]/1e3, 'k.', linewidth=2, markersize=8, label='VSMCa')

    plt.axvline(x=2.9, color='red', linestyle='--', linewidth=2)
    plt.axhline(y=16, color='red', linestyle='--', linewidth=2)
    plt.plot(2.9, 16, 'ro', linewidth=3, markersize=8)

    plt.title('Pressure vs Diameter')
    plt.xlabel('Diameter (mm)')
    plt.ylabel('Pressure (kPa)')
    plt.legend(loc='upper left')
    plt.ylim([0, 60])
    plt.grid(True)
    plt.show()



