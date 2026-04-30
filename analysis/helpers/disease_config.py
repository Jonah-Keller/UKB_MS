"""Loader for configs/disease.yaml — the disease-cohort swap configuration.

To replicate this study on a different cohort, edit configs/disease.yaml
(do not modify any analysis script). Every script that needs a disease-,
HLA-, or PRS-specific value imports the loaded config instead of hardcoding.

Usage:
    from analysis.helpers.disease_config import load_disease_config
    cfg = load_disease_config()
    icd = frozenset(cfg.icd_codes)
    cohort_short = cfg.cohort_short
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List

import yaml


def _project_root() -> Path:
    """Return the project root by walking up from this file."""
    return Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class DiseaseConfig:
    cohort_short: str
    disease_long: str
    disease_short_caps: str
    project_short: str
    project_title: str

    icd_codes: List[str]
    first_occurrence_field: str
    cohort_status_col: str
    status_values: Dict[str, str]
    control_exclusion_codes: Dict[str, List[str]]

    hla_allele: str
    hla_search_patterns: List[str]
    hla_carrier_col: str
    hla_dosage_col: str

    prs_pgs_ids: List[str]
    prs_label: str
    prs_combined_col: str

    cns_injury_markers: List[str]
    top_presymptomatic_deps: List[str]

    comparison_cohort_short: str = ""

    @property
    def all_exclusion_codes(self) -> frozenset:
        out: List[str] = []
        for codes in self.control_exclusion_codes.values():
            out.extend(codes)
        return frozenset(out)


@lru_cache(maxsize=1)
def load_disease_config(path: str | Path | None = None) -> DiseaseConfig:
    """Load configs/disease.yaml. Cached after first call."""
    if path is None:
        path = _project_root() / "configs" / "disease.yaml"
    with open(path, "r") as fh:
        raw = yaml.safe_load(fh)
    return DiseaseConfig(**raw)


if __name__ == "__main__":
    cfg = load_disease_config()
    print(f"cohort_short      = {cfg.cohort_short}")
    print(f"disease_long      = {cfg.disease_long}")
    print(f"icd_codes         = {cfg.icd_codes}")
    print(f"hla_allele        = {cfg.hla_allele}")
    print(f"prs_pgs_ids       = {cfg.prs_pgs_ids}")
    print(f"all_exclusion_codes ({len(cfg.all_exclusion_codes)}): "
          f"{sorted(cfg.all_exclusion_codes)[:5]}...")
