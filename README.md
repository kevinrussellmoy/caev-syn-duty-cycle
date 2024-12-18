# caev-syn-duty-cycle
Synthetic duty cycle generation from connected autonomous electric vehicle (C/AEV) driving data. This work supported the following publications, in reverse chronological order:

- Moy, K., Ganapathi, D., Geslin, A., Chueh, W., and Onori, S., “Synthetic duty cycles from real-world autonomous electric vehicle driving,” *Cell Reports Physical Science,* 2023, 101536.

General framework to run:
1. Run `syn_duty_cycle_gen.m` to generate all synthetic duty cycle information

-- Requires City 1 and City 2 drive cycle data as `City_1_CellPowerProfile.csv` and `City_2_CellPowerProfile.csv`. [The actual drive cycles are not supplied here as they are confidential. Any time-series data of velocity and battery current will do.]

-- User must select value for `city_select` to select City (`1`), City (`2`), or (`3`) City 1 and 2 combined.

-- User will obtain duty cycle current (in C-rate, normalized to nominal cell capacity `Q_nom`) and velocity.

-- Helper Functions:

---- `pca_k_means.m` to run the PCA + k-means synthetic duty cycle algorithm

-------- `mean_centering.m`, `start_end_disp.m`, and `rest_lengths.m` are helper functions for `pca_k_means.m`

2. Run `syn_duty_cycle_format` to plot the synthetic duty cycles and format/save them for later use in experimental test protocols.

-- Requires City 1 and City 2 drive cycle data as `City_1_CellPowerProfile.csv` and `City_2_CellPowerProfile.csv`. [The actual drive cycles are not supplied here as they are confidential. Any time-series data of velocity and battery current will do.]

3. Run `comp_ECAV_EV` for comparison between C/AEV driving and real-world electric vehicle driving

-- Requires City 1 and City 2 drive cycle data as `City_1_CellPowerProfile.csv` and `City_2_CellPowerProfile.csv`. [The actual drive cycles are not supplied here as they are confidential. Any time-series data of velocity and battery current will do.]

-- Requires downloading Vehicle Energy Dataset Dynamic Data (under \Data in [https://github.com/gsoh/VED/](https://github.com/gsoh/VED))
