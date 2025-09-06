import streamlit as st
import pandas as pd
from chembl_webresource_client.new_client import new_client
import plotly.express as px
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import time

# --- Constants ---
MAX_DRUGS = 75
TARGET_LIMIT = 3
WORKERS = 2
CACHE_SIZE = 2

# --- Session State ---
if 'results' not in st.session_state:
    st.session_state.results = {
        'drugs': pd.DataFrame(columns=['Drug', 'IC50', 'Phase', 'Target']),
        'gene': '',
        'threshold': 0,
        'last_updated': 0,
        'api_errors': 0
    }

# --- Debugged Target Fetching ---
@lru_cache(maxsize=CACHE_SIZE)
def get_targets(gene):
    try:
        # Try both synonym and pref_name searches with proper field filtering
        targets = list(new_client.target.filter(
            target_synonym__iexact=gene.upper()
        ).only('target_chembl_id', 'pref_name')[:TARGET_LIMIT*2])
        
        if not targets:
            targets = list(new_client.target.filter(
                pref_name__iexact=gene.upper()
            ).only('target_chembl_id', 'pref_name')[:TARGET_LIMIT*2])
        
        return targets if targets else None
    
    except Exception as e:
        st.session_state.results['api_errors'] += 1
        st.error(f"Target API error: {str(e)}")
        return None

# --- Debugged Drug Fetching ---
def fetch_drugs(target, threshold):
    drugs = []
    try:
        # Get mechanisms with proper field selection
        mechs = list(new_client.mechanism.filter(
            target_chembl_id=target['target_chembl_id']
        ).only('molecule_chembl_id', 'mechanism_of_action')[:MAX_DRUGS*2])
        
        for mech in mechs[:MAX_DRUGS//TARGET_LIMIT]:
            try:
                # Get minimal drug info
                drug = new_client.molecule.get(
                    mech['molecule_chembl_id'],
                    fields=['pref_name', 'max_phase', 'molecule_chembl_id']
                )
                
                # Get activities with proper filtering
                activities = list(new_client.activity.filter(
                    molecule_chembl_id=mech['molecule_chembl_id'],
                    target_chembl_id=target['target_chembl_id'],
                    standard_type="IC50",
                    standard_units="nM"
                ).only('standard_value', 'standard_units')[:1])
                
                if activities and activities[0]['standard_units'] == 'nM' and activities[0]['standard_value']:
                    ic50 = float(activities[0]['standard_value'])
                    if ic50 <= threshold:
                        drugs.append({
                            'Drug': drug.get('pref_name', f"Unknown_{mech['molecule_chembl_id'][:4]}"),
                            'IC50': ic50,
                            'Phase': drug.get('max_phase', 'Unknown'),
                            'Target': target.get('pref_name', target['target_chembl_id'])
                        })
            
            except Exception as e:
                st.warning(f"Skipped drug {mech.get('molecule_chembl_id')}: {str(e)}")
                continue
    
    except Exception as e:
        st.error(f"Mechanism API error for {target.get('target_chembl_id')}: {str(e)}")
    
    return drugs

# --- Streamlit UI ---
st.set_page_config(layout="wide")
st.title("💊 Debugged Drug Explorer")

col1, col2 = st.columns(2)
with col1:
    gene = st.text_input("Gene Symbol", "PTGS2")
with col2:
    threshold = st.number_input("IC50 Threshold (nM)", 
                              min_value=0.1, 
                              value=10000.0,
                              step=100.0)

if st.button("🔍 Analyze"):
    with st.spinner(f"Searching {gene} with threshold ≤ {threshold} nM..."):
        st.session_state.results['api_errors'] = 0  # Reset error counter
        
        # Get targets with debug info
        targets = get_targets(gene)
        if not targets:
            st.error(f"No targets found for {gene}. Try alternative gene names.")
            st.stop()
        
        st.info(f"Found {len(targets)} targets. Processing...")
        
        # Parallel processing with progress
        results = []
        with ThreadPoolExecutor(max_workers=WORKERS) as executor:
            futures = {executor.submit(fetch_drugs, t, threshold): t for t in targets[:TARGET_LIMIT]}
            for i, future in enumerate(as_completed(futures)):
                results.extend(future.result())
                st.write(f"Processed target {i+1}/{min(len(targets), TARGET_LIMIT)} | Found {len(results)} drugs")
                if len(results) >= MAX_DRUGS:
                    break
        
        # Handle results
        if not results:
            st.warning(f"No drugs found for {gene} at ≤ {threshold} nM. Try:")
            st.markdown("- Increasing threshold (current: {threshold} nM)")
            st.markdown("- Checking gene spelling (used: {gene.upper()})")
            st.markdown("- Reviewing API status: [ChEMBL Status](https://www.ebi.ac.uk/chembl/)")
            st.stop()
        
        df = pd.DataFrame(results).sort_values('IC50')
        st.session_state.results = {
            'drugs': df,
            'gene': gene,
            'threshold': threshold,
            'last_updated': time.time()
        }

# --- Display Results ---
if not st.session_state.results['drugs'].empty:
    st.success(f"Results for {st.session_state.results['gene']} (≤{st.session_state.results['threshold']} nM)")
    
    # Debug preview
    with st.expander("Raw Data Preview"):
        st.write(st.session_state.results['drugs'])
    
    # Safe visualization
    try:
        fig = px.scatter(
            st.session_state.results['drugs'],
            x='Drug',
            y='IC50',
            color='Phase',
            log_y=True,
            hover_data=['Target'],
            title=f"{st.session_state.results['gene']} Inhibitors"
        )
        st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.error(f"Visualization error: {str(e)}")
        st.write("Debug data:", st.session_state.results['drugs'])
    
    # Download
    st.download_button(
        "💾 Download CSV",
        st.session_state.results['drugs'].to_csv(index=False),
        f"{st.session_state.results['gene']}_drugs.csv"
    )

# --- System Status ---
st.caption(f"API errors in this session: {st.session_state.results.get('api_errors', 0)}")
if st.session_state.results.get('api_errors', 0) > 3:
    st.error("High API error rate. Consider pausing requests.")

# --- Cache Management ---
if st.session_state.results.get('last_updated', 0) < time.time() - 300:
    get_targets.cache_clear()
