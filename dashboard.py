import streamlit as st

from aneurysm_sim.model.model import simulate_aneurysm
from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import plots

@st.cache_resource
def run_simulation(genotype, treatment):
    params = ArterialParameters()
    return simulate_aneurysm(params, genotype=genotype, treatment=treatment)

def main():
    st.set_page_config(layout="wide", page_title="Aneurysm Simulation Dashboard")

    # Sidebar
    with st.sidebar:
        st.header("Aneurysm Simulation Dashboard")
        genotype = st.selectbox("Genotype", ["TT", "TC", "CC"])
        treatment = st.checkbox("Apply TGF-Beta Treatment", value=False)

    # Run main simulation
    results = run_simulation(genotype, treatment)
    # risk_score = calculate_risk_score(results)

    # Layout: plots side by side
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Normalised Densities Over Time")
        plt1 = plots.plot_normalised_densities(results)
        st.pyplot(plt1)

    with col2:
        st.subheader("Systolic Stretch Over Time")
        results_tt = run_simulation("TT", treatment)
        results_tc = run_simulation("TC", treatment)
        results_cc = run_simulation("CC", treatment)
        plt2 = plots.plot_systolic_stretch_over_time(results_tt, results_tc, results_cc)
        st.pyplot(plt2)

    # # Risk Score
    # st.markdown("---")
    # st.subheader("Risk Evaluation")
    # gauge_col1, gauge_col2, gauge_col3 = st.columns([1, 2, 1])
    # with gauge_col2:
    #     gauge = circular_gauge(risk_score)
    #     st.pyplot(gauge)

if __name__ == "__main__":
    main()