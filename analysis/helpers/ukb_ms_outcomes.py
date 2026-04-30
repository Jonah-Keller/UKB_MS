"""DEPRECATED back-compat shim: re-exports from ukb_disease_outcomes.

This module name is retained so that downstream stages refactored to read
from configs/disease.yaml in a later pass continue to work in the meantime.
New code should import from analysis.helpers.ukb_disease_outcomes directly.

The MS-specific names below are aliases of the generic disease-config-driven
names — they reflect whatever cohort is configured in configs/disease.yaml.
"""
from analysis.helpers.ukb_disease_outcomes import (
    DISEASE_ICD as MS_ICD,
    DISEASE_FIRST_OCCURRENCE_FIELD as MS_FIRST_OCCURRENCE_FIELD,
    DISEASE_EXCLUDE_CONTROL_CODES as MS_EXCLUDE_CONTROL_CODES,
    DEMYELINATING_CODES,
    EPILEPSY_CODES,
    MOVEMENT_CODES,
    DEMENTIA_CODES,
    is_disease_case as is_ms_case,
    is_valid_disease_control as is_valid_ms_control,
    classify_disease_timing as classify_ms_timing,
    years_to_disease_diagnosis as years_to_ms_diagnosis,
    CNS_INJURY_MARKERS,
    TOP_PRESYMPTOMATIC_DEPS,
)

__all__ = [
    "MS_ICD",
    "MS_FIRST_OCCURRENCE_FIELD",
    "MS_EXCLUDE_CONTROL_CODES",
    "DEMYELINATING_CODES",
    "EPILEPSY_CODES",
    "MOVEMENT_CODES",
    "DEMENTIA_CODES",
    "is_ms_case",
    "is_valid_ms_control",
    "classify_ms_timing",
    "years_to_ms_diagnosis",
    "CNS_INJURY_MARKERS",
    "TOP_PRESYMPTOMATIC_DEPS",
]
