import streamlit as st

from aneurysm_sim.model.model import simulate_aneurysm
from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import plots

def main():
    st.title("Aneurysm Simulation Dashboard")
    
    # Sidebar for user inputs
    st.sidebar.header("Simulation Parameters")
    genotype = st.sidebar.selectbox("Select Genotype", ["TT", "TC", "CC"])
    treatment = st.sidebar.checkbox("Apply TGF-Beta Treatment", value=False)
    
    # Run simulation
    params = ArterialParameters()
    results = simulate_aneurysm(params, genotype=genotype, treatment=treatment)
    
    # Display plots
    st.subheader("Normalised Densities Over Time")
    plt1 = plots.plot_normalised_densities(results)
    st.pyplot(plt1)
    
    st.subheader("Systolic Stretch Over Time by Genotype")
    results_tt = simulate_aneurysm(params, genotype="TT", treatment=treatment)
    results_tc = simulate_aneurysm(params, genotype="TC", treatment=treatment)
    results_cc = simulate_aneurysm(params, genotype="CC", treatment=treatment)
    plt2 = plots.plot_systolic_stretch_over_time(results_tt, results_tc, results_cc)
    st.pyplot(plt2)

if __name__ == "__main__":
    main()