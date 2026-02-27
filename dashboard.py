'''This is a work in progress Streamlit dashboard for aneurysm simulation.'''

import streamlit as st
import plotly.graph_objects as go

from aneurysm_sim.model.model import simulate_aneurysm
from aneurysm_sim.config.parameters import ArterialParameters
from aneurysm_sim.model import plots

def circular_gauge(abr_score):
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=abr_score,
        domain={'x': [0, 1], 'y': [0, 1]},
        title={'text': "ABR (Wall Stress / Strength)", 'font': {'size': 14}},  # smaller title
        gauge={
            'axis': {'range': [0, 2], 'tickwidth': 1, 'tickcolor': "darkgray"},
            'bar': {'color': "darkblue"},
            'steps': [
                {'range': [0, 0.4], 'color': "lightgreen"},
                {'range': [0.41, 0.7], 'color': "gold"},
                {'range': [0.71, 10.0], 'color': "red"},
            ],
            'threshold': {
                'line': {'color': "black", 'width': 4},
                'thickness': 0.75,
                'value': 1.0
            }
        }
    ))
    fig.update_layout(
        margin=dict(t=10, b=10, l=10, r=10),
        width=300,   # set smaller width
        height=250   # set smaller height
    )
    return fig

@st.cache_resource
def run_simulation(genotype, polygenic_score, treatment):
    params = ArterialParameters(genotype=genotype, polygenic_score=polygenic_score)
    return simulate_aneurysm(params, treatment=treatment)

def main():
    st.set_page_config(layout="wide", page_title="Aneurysm Simulation Dashboard")

    # Sidebar
    with st.sidebar:
        st.header("Aneurysm Simulation Dashboard")
        genotype = st.selectbox("Genotype", ["TT", "TC", "CC"])
        polygenic_score = st.selectbox("Polygenic Score", [0, 1, 2, 3, 4])
        treatment = st.checkbox("Apply TGF-Beta Treatment", value=False)

    # Run main simulation
    results = run_simulation(genotype, polygenic_score, treatment)
    # risk_score = calculate_risk_score(results)

    # Layout: plots side by side
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Normalised Densities Over Time")
        fig1 = plots.plot_normalised_densities(results, legend=True)
        st.pyplot(fig1)

    with col2:
        st.subheader("Systolic Stretch Over Time")
        results_tt = run_simulation("TT", polygenic_score, treatment)
        results_tc = run_simulation("TC", polygenic_score, treatment)
        results_cc = run_simulation("CC", polygenic_score, treatment)
        fig2 = plots.plot_stretch_by_genotype(results_tt, results_tc, results_cc)
        st.pyplot(fig2)

    # # Risk Score
    # st.subheader("Aneurysm Risk Score (ABR)")
    # abr_score = results["abr_score"]  # Assuming abr_score is part of the results
    # st.plotly_chart(circular_gauge(abr_score), use_container_width=True)
    # st.write(f"**ABR Score:** {abr_score:.2f}")


if __name__ == "__main__":
    main()