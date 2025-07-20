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
        title={'text': "ABR (Wall Stress / Strength)", 'font': {'size': 18}},
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
    fig.update_layout(margin=dict(t=0, b=0, l=0, r=0))
    return fig

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

    # Risk Score
    st.subheader("Aneurysm Risk Score (ABR)")
    abr_score = results["abr_score"]  # Assuming abr_score is part of the results
    st.plotly_chart(circular_gauge(abr_score), use_container_width=True)
    st.write(f"**ABR Score:** {abr_score:.2f}")


if __name__ == "__main__":
    main()