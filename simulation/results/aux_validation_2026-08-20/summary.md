# Aux-pipeline validation -- 2026-08-20

12/29 gated checks passed (16 rows tracked/informational, not gated).

| scenario | category | status | pass | err_cells | notes |
|---|---|---|---|---|---|
| baseline_measured_conditions | position | ok | False | -17.893106792034892 |  |
| baseline_measured_conditions | position | ok | False | 16.436124204544758 |  |
| baseline_measured_conditions | position | ok | False | -18.473487733078493 |  |
| baseline_measured_conditions | position | ok | False | -0.09239339827795787 |  |
| baseline_tool_default | position | ok | True | 0.35139089734882617 |  |
| baseline_tool_default | position | ok | True | -0.23005679190908643 |  |
| baseline_tool_default | position | ok | True | 0.19480020790629426 |  |
| baseline_tool_default | position | ok | True | -0.20037997329210444 |  |
| doublet_2cell | resolution | resolved | True |  | valley 8.7/5.4 dB below the two peaks |
| doublet_1cell | resolution | resolved | False |  | valley 25.3/15.8 dB below the two peaks |
| weak_near_floor | amplitude | ok | None | 0.35139089734882617 |  |
| weak_near_floor | amplitude | ok | None | -0.10018998664605222 |  |
| weak_near_floor | amplitude | ok | None | -0.15028497996516402 |  |
| weak_near_floor | amplitude | ok | None | -0.20037997329210444 |  |
| weak_near_floor_elevated_noise | amplitude | ok | None | 0.35139213716554835 |  |
| weak_near_floor_elevated_noise | amplitude | ok | None | -0.10016303411516207 |  |
| weak_near_floor_elevated_noise | amplitude | ok | None | -0.1502445511727431 |  |
| weak_near_floor_elevated_noise | amplitude | ok | None | -0.20032606823032414 |  |
| nyquist_edge_inside | range | ok | True | -0.48149747068475357 |  |
| nyquist_beyond_aliased | range | aliased_as_expected | True | 0.8047021187005436 |  |
| nyquist_beyond_aliased | range | aliased_as_expected | True | -1.5947075654345586 |  |
| ripple_stress_near_limit | robustness | ok | None | -0.7339411492554491 | 23.13% backward steps (expected at this ripple) |
| ripple_stress_near_limit | robustness | ok | None | 3.2568973407076753 | 23.13% backward steps (expected at this ripple) |
| ripple_stress_near_limit | robustness | ok | None | -17.10083022590454 | 23.13% backward steps (expected at this ripple) |
| ripple_stress_near_limit | robustness | ok | None | 0.04600438020792897 | 23.13% backward steps (expected at this ripple) |
| bow_stress_2x | robustness | ok | None | 0.2367225380126685 | 0.00% backward steps (expected at this ripple) |
| bow_stress_2x | robustness | ok | None | 0.022105061569630286 | 0.00% backward steps (expected at this ripple) |
| bow_stress_2x | robustness | ok | None | 1.577366184878791 | 0.00% backward steps (expected at this ripple) |
| bow_stress_2x | robustness | ok | None | -0.18596591606190252 | 0.00% backward steps (expected at this ripple) |
| tau_aux_miscal_0.5pct | calibration | ok | False | 16.488929697828446 | tau_pct=0.5, vs.predicted=+0.35 cells |
| tau_aux_miscal_0.5pct | calibration | ok | False | 105.01476147208517 | tau_pct=0.5, vs.predicted=-0.23 cells |
| tau_aux_miscal_0.5pct | calibration | ok | False | 368.55166413188294 | tau_pct=0.5, vs.predicted=+0.19 cells |
| tau_aux_miscal_0.5pct | calibration | ok | False | 701.4317417866965 | tau_pct=0.5, vs.predicted=-0.20 cells |
| tau_aux_miscal_1pct | calibration | ok | False | 32.46669088642454 | tau_pct=1.0, vs.predicted=+0.35 cells |
| tau_aux_miscal_1pct | calibration | ok | False | 209.21755183249482 | tau_pct=1.0, vs.predicted=-0.23 cells |
| tau_aux_miscal_1pct | calibration | ok | False | 733.2614303933233 | tau_pct=1.0, vs.predicted=+0.19 cells |
| tau_aux_miscal_1pct | calibration | ok | False | 1396.1170108560693 | tau_pct=1.0, vs.predicted=-0.20 cells |
| tau_aux_miscal_2pct | calibration | ok | False | 63.95227911100619 | tau_pct=2.0, vs.predicted=+0.35 cells |
| tau_aux_miscal_2pct | calibration | ok | False | 414.5583446015069 | tau_pct=2.0, vs.predicted=-0.23 cells |
| tau_aux_miscal_2pct | calibration | ok | False | 1451.954205084863 | tau_pct=2.0, vs.predicted=+0.19 cells |
| tau_aux_miscal_2pct | calibration | ok | False | 2765.0556293161517 | tau_pct=2.0, vs.predicted=-0.20 cells |
| dispersion_beta2_off | known_limitation | ok | True | -0.38017597230587846 |  |
| dispersion_beta2_off | known_limitation | ok | True | 1.0991170772456405 |  |
| dispersion_beta2_on | known_limitation | ok | True | -0.38017597230587846 |  |
| dispersion_beta2_on | known_limitation | ok | True | 1.0991170772456405 |  |
