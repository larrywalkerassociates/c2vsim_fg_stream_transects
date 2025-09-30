# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: editable,slideshow,tags,-all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
# ---

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Introduction
# - Problem Statement:
#    - 10 years of SGMA have generated large amounts of field and modeling data
#    - Presenting these data in a meaningful, intuitive, and transparent way requires sophisticated and highly customizable
#      visualization tools
#    - Mainstream tools are propietary, resulting in license costs and limited flexibility for users
#  - Why it is important:
#    - Field and simulated data is of limited use if it cannot be effectively communicated to stakeholders and used
#      for decision-making
#    - When processing data, the devil often lies in the details, which proprietary tools often obscure making it difficult
#      to peer-review and recreate workflows
#    - High-level programming languages (e.g., Python, R) provide a wide array of open-source and accessible tools that
#      can be used to streamline data processing, visualization, and dashboard deployment for engineers and scientists.
#  - How we have approached it (in general terms):
#    - We have developed a Python-based workflow to process public groundwater levels, lithology logs, bathymetry, and
#      simulated groundwater levels, displaying them in animated stream transects, which are deployed as dashboards.
#   - What is new about what we have done:
#     - Integrating multiple open-source datasets into one type of visualization
#     - Customizable and transparent workflow that can be replicated at other streams in the Central Valley with low
#       effort
