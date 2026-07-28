import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os
from sklearn.decomposition import PCA

# ---------------------------------------------------------
# 1. Page Configuration & Title
# ---------------------------------------------------------
st.set_page_config(page_title="Rice Microbiome Dashboard", layout="wide")
st.title("🌾 Rice Microbiome: Complete Ecological Analysis Dashboard")

# ---------------------------------------------------------
# 2. Resilient Data Loading & Processing Engine
# ---------------------------------------------------------
@st.cache_data
def load_and_clean_data():
    local_files = os.listdir('.')
    
    # Initialize variables
    pmoa_raw, mcra_raw, pmoa_major_raw, mcra_major_raw = None, None, None, None
    
    # 1. Try finding standalone CSV files first
    pmoa_rel_csv = next((f for f in local_files if 'pmoA_OTUs_rel_abundance' in f and f.endswith('.csv')), None)
    mcra_rel_csv = next((f for f in local_files if 'mcrA_OTUs_rel_abundance' in f and f.endswith('.csv')), None)
    pmoa_major_csv = next((f for f in local_files if 'Major_pmoA_OTUs' in f and f.endswith('.csv')), None)
    mcra_major_csv = next((f for f in local_files if 'Major_mcrA_OTUs' in f and f.endswith('.csv')), None)
    
    if pmoa_rel_csv and mcra_rel_csv and pmoa_major_csv and mcra_major_csv:
        pmoa_raw = pd.read_csv(pmoa_rel_csv)
        mcra_raw = pd.read_csv(mcra_rel_csv)
        pmoa_major_raw = pd.read_csv(pmoa_major_csv)
        mcra_major_raw = pd.read_csv(mcra_major_csv)
    else:
        # 2. Fallback to the Master Excel Workbook
        excel_file = next((f for f in local_files if 'IU_pmoA_mcrA_2025_Data' in f and f.endswith(('.xlsx', '.xls'))), None)
        if excel_file:
            try:
                pmoa_raw = pd.read_excel(excel_file, sheet_name="pmoA_OTUs_rel_abundance", engine="openpyxl")
                mcra_raw = pd.read_excel(excel_file, sheet_name="mcrA_OTUs_rel_abundance", engine="openpyxl")
                pmoa_major_raw = pd.read_excel(excel_file, sheet_name="Major_pmoA_OTUs", engine="openpyxl")
                mcra_major_raw = pd.read_excel(excel_file, sheet_name="Major_mcrA_OTUs", engine="openpyxl")
            except Exception as e:
                st.error(f"Error reading Excel: {e}")
                st.stop()
                
    if pmoa_raw is None or mcra_raw is None:
        st.error("### ❌ Could not find data files.")
        st.stop()
        
    # Column mapping logic to strictly enforce MAJOR OTUs only
    pmoa_id_col = 'OTU_No.' if 'OTU_No.' in pmoa_major_raw.columns else pmoa_major_raw.columns[0]
    mcra_id_col = 'OTU_No.' if 'OTU_No.' in mcra_major_raw.columns else mcra_major_raw.columns[0]
    
    pmoa_major_list = [f"Otu{int(otu.lower().split('otu')[1]):03d}" for otu in pmoa_major_raw[pmoa_id_col].astype(str).str.strip() if 'otu' in otu.lower()]
    mcra_major_list = [f"Otu{int(otu.lower().split('otu')[1]):02d}" for otu in mcra_major_raw[mcra_id_col].astype(str).str.strip() if 'otu' in otu.lower()]
            
    pmoa_major_cols = [col for col in pmoa_raw.columns if col in pmoa_major_list]
    mcra_major_cols = [col for col in mcra_raw.columns if col in mcra_major_list]
    
    pmoa_df = pmoa_raw[['Group', 'Stage'] + pmoa_major_cols].copy()
    mcra_df = mcra_raw[['Group', 'Stage'] + mcra_major_cols].copy()
    
    # Clean and parse metadata
    def parse_metadata(df):
        varieties, treatments = [], []
        for group in df['Group']:
            if 'Before_flooding' in str(group) or 'Before flooding' in str(group):
                varieties.append('Both')
                treatments.append('Baseline')
            else:
                parts = str(group).split('_')
                varieties.append(parts[0])
                treatments.append(parts[1])
        df['Variety'] = varieties
        df['Treatment'] = treatments
        return df

    pmoa_df = parse_metadata(pmoa_df)
    mcra_df = parse_metadata(mcra_df)
    
    # Establish chronological order
    stage_order = {'Before flooding': 1, 'Early_tillering_stage': 2, 'Panicle_formation_stage': 3, 'Early_heading_stage': 4}
    pmoa_df['Stage_Order'] = pmoa_df['Stage'].map(stage_order)
    mcra_df['Stage_Order'] = mcra_df['Stage'].map(stage_order)
    
    return pmoa_df, mcra_df, pmoa_major_cols, mcra_major_cols

def calculate_shannon(row, cols):
    vals = row[cols].astype(float).values
    vals = vals[vals > 0]
    if len(vals) == 0: return 0
    p = vals / vals.sum()
    return -np.sum(p * np.log(p))

# Load the data
pmoa_df, mcra_df, pmoa_cols, mcra_cols = load_and_clean_data()

# Calculate Diversity and Ratios
pmoa_df['pmoA_Shannon'] = pmoa_df.apply(lambda r: calculate_shannon(r, pmoa_cols), axis=1)
mcra_df['mcrA_Shannon'] = mcra_df.apply(lambda r: calculate_shannon(r, mcra_cols), axis=1)

merged_df = pd.merge(
    pmoa_df[['Group', 'Stage', 'Variety', 'Treatment', 'Stage_Order', 'pmoA_Shannon']],
    mcra_df[['Group', 'mcrA_Shannon']], on='Group'
)
merged_df['Tug_of_War_Ratio'] = np.log2(merged_df['pmoA_Shannon'] / merged_df['mcrA_Shannon'])

# ---------------------------------------------------------
# 3. Dashboard Tabs Definition
# ---------------------------------------------------------
tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs([
    "Graph A: Tug-of-War", 
    "Graph B: Succession", 
    "Graph C: Diversity", 
    "Graph D: Microbe Inspector", 
    "Graph E & Stats", 
    "Community Snapshot",
    "Graph F & G: Shift & Radar", 
    "Graph H & I: Advanced Ecology"
])

# Global styling and order rules
color_map = {"Control": "#ef553b", "KH32C": "#00cc96", "Baseline": "#636efa"}
stage_cats = ["Before flooding", "Early_tillering_stage", "Panicle_formation_stage", "Early_heading_stage"]

# ---------------------------------------------------------
# Tab 1: Core Methane Tug-of-War
# ---------------------------------------------------------
with tab1:
    st.header("Graph A: Core Methane Tug-of-War")
    baseline_data = merged_df[merged_df['Treatment'] == 'Baseline']
    treated_data = merged_df[merged_df['Treatment'] != 'Baseline']
    recon_list = []
    
    for var in ['IR64', 'Nipponbare']:
        for trt in ['Control', 'KH32C']:
            b_copy = baseline_data.copy()
            b_copy['Variety'] = var
            b_copy['Treatment'] = trt
            recon_list.append(b_copy)
    recon_list.append(treated_data)
    
    if recon_list:
        final_line_df = pd.concat(recon_list).groupby(['Stage', 'Variety', 'Treatment', 'Stage_Order']).mean(numeric_only=True).reset_index().sort_values('Stage_Order')
        fig_a = px.line(final_line_df, x='Stage', y='Tug_of_War_Ratio', color='Treatment', facet_col='Variety', markers=True, color_discrete_map=color_map, category_orders={"Stage": stage_cats}, title="Log2 Ratio of Methane Eaters vs Producers")
        st.plotly_chart(fig_a, use_container_width=True, key="tug_of_war_chart")

# ---------------------------------------------------------
# Tab 2: Top-Driver Community Succession
# ---------------------------------------------------------
with tab2:
    st.header("Graph B: Top-Driver Community Succession")
    
    def plot_stacked_bars(df, cols, title, chart_key):
        melted = df.melt(id_vars=['Group', 'Stage', 'Variety', 'Treatment', 'Stage_Order'], value_vars=cols, var_name='OTU', value_name='Abundance')
        avg_melted = melted.groupby(['Stage', 'Variety', 'Treatment', 'Stage_Order', 'OTU']).mean(numeric_only=True).reset_index()
        top_5 = avg_melted.groupby('OTU')['Abundance'].sum().nlargest(5).index
        avg_melted['OTU_Group'] = avg_melted['OTU'].apply(lambda x: x if x in top_5 else 'Others')
        final_bar = avg_melted.groupby(['Stage', 'Variety', 'Treatment', 'Stage_Order', 'OTU_Group']).sum(numeric_only=True).reset_index().sort_values('Stage_Order')
        
        if not final_bar.empty:
            fig = px.bar(final_bar[final_bar['Treatment'] != 'Baseline'], x='Stage', y='Abundance', color='OTU_Group', facet_col='Variety', facet_row='Treatment', title=title, category_orders={"Stage": stage_cats})
            st.plotly_chart(fig, use_container_width=True, key=chart_key)
        
    plot_stacked_bars(pmoa_df, pmoa_cols, "pmoA (Methane Eaters)", "pmoa_bar_chart")
    plot_stacked_bars(mcra_df, mcra_cols, "mcrA (Methane Producers)", "mcra_bar_chart")

# ---------------------------------------------------------
# Tab 3: Microbial Alpha Diversity
# ---------------------------------------------------------
with tab3:
    st.header("Graph C: Microbial Alpha Diversity Tracking")
    c1, c2 = st.columns(2)
    with c1:
        fig_box_p = px.box(pmoa_df[pmoa_df['Treatment'] != 'Baseline'], x='Stage', y='pmoA_Shannon', color='Treatment', facet_col='Variety', points="all", title="pmoA Alpha Diversity", color_discrete_map=color_map, category_orders={"Stage": stage_cats})
        st.plotly_chart(fig_box_p, use_container_width=True, key="alpha_pmoa_box")
    with c2:
        fig_box_m = px.box(mcra_df[mcra_df['Treatment'] != 'Baseline'], x='Stage', y='mcrA_Shannon', color='Treatment', facet_col='Variety', points="all", title="mcrA Alpha Diversity", color_discrete_map=color_map, category_orders={"Stage": stage_cats})
        st.plotly_chart(fig_box_m, use_container_width=True, key="alpha_mcra_box")

# ---------------------------------------------------------
# Tab 4: Individual Core Microbe Tracker
# ---------------------------------------------------------
with tab4:
    st.header("Graph D: Individual Core Microbe Tracker")
    c1, c2 = st.columns(2)
    with c1: selected_pmoa = st.selectbox("Choose Major pmoA:", options=pmoa_cols, key="s1")
    with c2: selected_mcra = st.selectbox("Choose Major mcrA:", options=mcra_cols, key="s2")
    
    def fill_baseline(df, col_name):
        base = df[df['Treatment'] == 'Baseline']
        treat = df[df['Treatment'] != 'Baseline']
        blist = []
        for v in ['IR64', 'Nipponbare']:
            for t in ['Control', 'KH32C']:
                bc = base.copy()
                bc['Variety'] = v
                bc['Treatment'] = t
                blist.append(bc)
        blist.append(treat)
        return pd.concat(blist).groupby(['Stage', 'Variety', 'Treatment', 'Stage_Order']).mean(numeric_only=True).reset_index().sort_values('Stage_Order')

    st.plotly_chart(px.line(fill_baseline(pmoa_df, selected_pmoa), x='Stage', y=selected_pmoa, color='Treatment', facet_col='Variety', markers=True, color_discrete_map=color_map, category_orders={"Stage": stage_cats}), use_container_width=True, key="line1")
    st.plotly_chart(px.line(fill_baseline(mcra_df, selected_mcra), x='Stage', y=selected_mcra, color='Treatment', facet_col='Variety', markers=True, color_discrete_map=color_map, category_orders={"Stage": stage_cats}), use_container_width=True, key="line2")

# ---------------------------------------------------------
# Tab 5: Heatmap & Beta Diversity (PCA)
# ---------------------------------------------------------
with tab5:
    st.header("Graph E & Stats: Heatmap & Beta Diversity")
    threshold = st.slider("Min Abundance Threshold (%) for Heatmap", 0.0, 10.0, 1.0, step=0.5) / 100.0
    
    pmoa_subset = pmoa_df[pmoa_df['Treatment'] != 'Baseline']
    if not pmoa_subset.empty:
        pmoa_avg = pmoa_subset.groupby(['Variety', 'Treatment', 'Stage'])[pmoa_cols].mean().reset_index()
        pmoa_melt = pmoa_avg.melt(id_vars=['Variety', 'Treatment', 'Stage'], var_name='OTU', value_name='Abundance')
        core_otus = pmoa_melt.groupby('OTU')['Abundance'].max()
        valid_otus = core_otus[core_otus >= threshold].index
        
        if len(valid_otus) > 0:
            fig_heat = px.density_heatmap(pmoa_melt[pmoa_melt['OTU'].isin(valid_otus)], x='Stage', y='OTU', z='Abundance', facet_col='Treatment', facet_row='Variety', color_continuous_scale='Viridis', title="pmoA Core Heatmap")
            st.plotly_chart(fig_heat, use_container_width=True, key="heatmap_pmoa")
    
    st.subheader("Beta Diversity (PCA Mapping)")
    pca_data = pmoa_df[pmoa_df['Treatment'] != 'Baseline'].copy()
    if len(pca_data) > 3: 
        pca = PCA(n_components=2)
        components = pca.fit_transform(pca_data[pmoa_cols].fillna(0))
        pca_data['PCA1'], pca_data['PCA2'] = components[:, 0], components[:, 1]
        fig_pca = px.scatter(pca_data, x='PCA1', y='PCA2', color='Treatment', symbol='Stage', facet_col='Variety', title="PCA Beta Diversity", color_discrete_map=color_map)
        st.plotly_chart(fig_pca, use_container_width=True, key="pca_beta_div")

# ---------------------------------------------------------
# Tab 6: Community Snapshot
# ---------------------------------------------------------
with tab6:
    st.header("Community Snapshot")
    chosen_stage = st.selectbox("Select Development Stage:", options=stage_cats, key="snapshot_stage")
    
    if chosen_stage == 'Before flooding':
        stage_pmoa = pmoa_df[pmoa_df['Stage'] == chosen_stage]
        stage_mcra = mcra_df[mcra_df['Stage'] == chosen_stage]
    else:
        chosen_variety = st.radio("Select Variety:", options=['IR64', 'Nipponbare'], key="snapshot_variety")
        stage_pmoa = pmoa_df[(pmoa_df['Stage'] == chosen_stage) & (pmoa_df['Variety'] == chosen_variety)]
        stage_mcra = mcra_df[(mcra_df['Stage'] == chosen_stage) & (mcra_df['Variety'] == chosen_variety)]
    
    if stage_pmoa.empty or stage_mcra.empty:
        st.warning(f"No data available for {chosen_stage}")
    else:
        pmoa_melt = stage_pmoa.melt(id_vars=['Treatment'], value_vars=pmoa_cols, var_name='OTU', value_name='Relative Abundance')
        mcra_melt = stage_mcra.melt(id_vars=['Treatment'], value_vars=mcra_cols, var_name='OTU', value_name='Relative Abundance')

        c1, c2 = st.columns(2)
        with c1: st.plotly_chart(px.bar(pmoa_melt.groupby(['Treatment', 'OTU']).mean().reset_index(), x='OTU', y='Relative Abundance', color='Treatment', barmode='group', title="pmoA Profile"), use_container_width=True, key="snap_p")
        with c2: st.plotly_chart(px.bar(mcra_melt.groupby(['Treatment', 'OTU']).mean().reset_index(), x='OTU', y='Relative Abundance', color='Treatment', barmode='group', title="mcrA Profile"), use_container_width=True, key="snap_m")

# ---------------------------------------------------------
# Tab 7: THE BIOFERTILIZER SHIFT (FOLD CHANGE & RADAR)
# ---------------------------------------------------------
with tab7:
    st.header("Graph F & G: The Biofertilizer Shift")
    st.markdown("Calculates the exact difference between the Control and KH32C mud to explicitly prove how the community shifted.")
    
    shift_stage = st.selectbox("Select Stage to Analyze Shift:", options=["Early_tillering_stage", "Panicle_formation_stage", "Early_heading_stage"], key="shift_stage")
    shift_var = st.radio("Select Variety to Analyze Shift:", options=["IR64", "Nipponbare"], key="shift_var_radio")
    
    st.markdown("---")
    
    def plot_fold_change(df, cols, title_prefix, chart_key):
        stage_df = df[(df['Stage'] == shift_stage) & (df['Variety'] == shift_var)]
        if stage_df.empty: return
        
        eps = 1e-5 
        ctrl_mean = stage_df[stage_df['Treatment'] == 'Control'][cols].mean(numeric_only=True) + eps
        kh_mean = stage_df[stage_df['Treatment'] == 'KH32C'][cols].mean(numeric_only=True) + eps
        
        l2fc = np.log2(kh_mean / ctrl_mean).reset_index()
        l2fc.columns = ['OTU', 'Log2_Fold_Change']
        
        l2fc['Shift_Direction'] = l2fc['Log2_Fold_Change'].apply(lambda x: 'Enriched by KH32C (Green)' if x > 0 else 'Suppressed by KH32C (Red)')
        l2fc = l2fc.sort_values('Log2_Fold_Change')
        
        fig = px.bar(l2fc, x='Log2_Fold_Change', y='OTU', orientation='h', color='Shift_Direction',
                     title=f"Graph F: {title_prefix} Fold-Change",
                     color_discrete_map={'Enriched by KH32C (Green)': '#00cc96', 'Suppressed by KH32C (Red)': '#ef553b'})
        fig.add_vline(x=0, line_width=2, line_dash="dash", line_color="black")
        st.plotly_chart(fig, use_container_width=True, key=chart_key)

    c1, c2 = st.columns(2)
    with c1: plot_fold_change(pmoa_df, pmoa_cols, "Methane Eaters", "fc_pmoa")
    with c2: plot_fold_change(mcra_df, mcra_cols, "Methane Producers", "fc_mcra")

    st.markdown("---")
    st.subheader("Graph G: Microbial Fingerprint (Radar Chart)")
    
    def plot_radar(df, cols, title_prefix, chart_key):
        stage_df = df[(df['Stage'] == shift_stage) & (df['Variety'] == shift_var)]
        if stage_df.empty: return
        
        top_8 = stage_df[cols].mean().nlargest(8).index.tolist()
        ctrl_vals = stage_df[stage_df['Treatment'] == 'Control'][top_8].mean().tolist()
        kh_vals = stage_df[stage_df['Treatment'] == 'KH32C'][top_8].mean().tolist()
        
        top_8.append(top_8[0])
        ctrl_vals.append(ctrl_vals[0])
        kh_vals.append(kh_vals[0])
        
        fig = go.Figure()
        fig.add_trace(go.Scatterpolar(r=ctrl_vals, theta=top_8, fill='toself', name='Control', line_color='#ef553b'))
        fig.add_trace(go.Scatterpolar(r=kh_vals, theta=top_8, fill='toself', name='KH32C', line_color='#00cc96'))
        
        fig.update_layout(polar=dict(radialaxis=dict(visible=True)), showlegend=True, title=f"{title_prefix} Community Shape")
        st.plotly_chart(fig, use_container_width=True, key=chart_key)

    c3, c4 = st.columns(2)
    with c3: plot_radar(pmoa_df, pmoa_cols, "pmoA", "radar_pmoa")
    with c4: plot_radar(mcra_df, mcra_cols, "mcrA", "radar_mcra")

# ---------------------------------------------------------
# Tab 8: ADVANCED ECOLOGICAL DYNAMICS
# ---------------------------------------------------------
with tab8:
    st.header("Graph H: The Microbial River (Continuous Flow)")
    st.markdown("Visualizes the community as a fluid ecosystem. Watch for 'blooms' (sudden expansions) in the stream across growth stages.")
    
    river_var = st.radio("Select Variety for the River:", options=["IR64", "Nipponbare"], key="river_var", horizontal=True)
    
    # 1. Create a reusable function to build the River Graph
    def plot_microbial_river(df, cols, title_prefix, chart_key):
        river_data = df[(df['Variety'] == river_var) & (df['Treatment'] != 'Baseline')]
        
        if not river_data.empty:
            river_melt = river_data.melt(id_vars=['Stage', 'Stage_Order', 'Treatment'], value_vars=cols, var_name='OTU', value_name='Abundance')
            river_avg = river_melt.groupby(['Stage', 'Stage_Order', 'Treatment', 'OTU']).mean().reset_index().sort_values('Stage_Order')
            
            # Keep top 8 for a smooth visual, group the rest
            top_river = river_avg.groupby('OTU')['Abundance'].sum().nlargest(8).index
            river_avg['OTU_Group'] = river_avg['OTU'].apply(lambda x: x if x in top_river else f'Other_{title_prefix}')
            final_river = river_avg.groupby(['Stage', 'Stage_Order', 'Treatment', 'OTU_Group']).sum(numeric_only=True).reset_index().sort_values('Stage_Order')
            
            fig_river = px.area(final_river, x='Stage', y='Abundance', color='OTU_Group', facet_col='Treatment', 
                                line_shape='spline', # Creates the fluid river effect
                                title=f"{title_prefix} Microbial River Flow ({river_var})",
                                category_orders={"Stage": stage_cats})
            st.plotly_chart(fig_river, use_container_width=True, key=chart_key)
        else:
            st.warning(f"Not enough data to map the {title_prefix} microbial river.")

    # 2. Render both rivers sequentially for full-width visibility
    plot_microbial_river(pmoa_df, pmoa_cols, "pmoA (Methane Eaters)", "river_pmoa")
    plot_microbial_river(mcra_df, mcra_cols, "mcrA (Methane Producers)", "river_mcra")

    st.markdown("---")
    
    st.header("Graph I: Ecological Cross-Talk (The Microbe Matrix)")
    st.markdown("**Blue = Cooperation** (growing together), **Red = Competition** (fighting for resources).")
    
    matrix_trt = st.selectbox("Select Treatment to Analyze Cross-Talk:", options=["Control", "KH32C"], key="matrix_trt")
    matrix_var = st.selectbox("Select Variety for Cross-Talk:", options=["IR64", "Nipponbare"], key="matrix_var_select")
    
    p_sub = pmoa_df[(pmoa_df['Treatment'] == matrix_trt) & (pmoa_df['Variety'] == matrix_var)]
    m_sub = mcra_df[(mcra_df['Treatment'] == matrix_trt) & (mcra_df['Variety'] == matrix_var)]
    
    if len(p_sub) > 2 and len(m_sub) > 2:
        top_p = p_sub[pmoa_cols].mean().nlargest(5).index.tolist()
        top_m = m_sub[mcra_cols].mean().nlargest(5).index.tolist()
        
        cross_df = pd.merge(p_sub[['Group'] + top_p], m_sub[['Group'] + top_m], on='Group')
        
        corr_matrix = cross_df[top_p + top_m].corr().loc[top_p, top_m]
        
        fig_corr = px.imshow(corr_matrix, 
                             labels=dict(x="mcrA (Methane Producers)", y="pmoA (Methane Eaters)", color="Correlation"),
                             x=top_m, y=top_p,
                             color_continuous_scale='RdBu', zmin=-1, zmax=1,
                             title=f"Cross-Talk Matrix: {matrix_trt} on {matrix_var}")
        
        fig_corr.update_traces(text=corr_matrix.round(2).values, texttemplate="%{text}")
        st.plotly_chart(fig_corr, use_container_width=True, key="cross_talk_matrix")
    else:
        st.warning("Need more sample replicates to calculate statistical correlation.")