#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd

# ============================================================
# LAB Genome Landscape
# Reproducible downstream analysis of fermentation markers
# ============================================================

ROOT = Path(__file__).resolve().parents[2]

METADATA = ROOT / "metadata/processed/genome_metadata.csv"
UNIQUE_MARKERS = ROOT / "metadata/processed/fermentation_unique_markers.csv"
GENE_MATRIX = ROOT / "analysis/functional_modules/fermentation/fermentation_gene_matrix.csv"

OUTDIR = ROOT / "results/reproducibility/recomputed_v1/tables"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Operational thresholds used in the manuscript
CORE_THRESHOLD = 1.00
NEAR_CORE_THRESHOLD = 0.95
VARIABLE_THRESHOLD = 0.15

MODULE_NAMES = {
    "FM01": "EMP_glycolysis",
    "FM02": "Lactate_formation",
    "FM03": "Phosphoketolase_pathway",
    "FM04": "Acetate_branch",
    "FM05": "Ethanol_branch",
    "FM06": "Pyruvate_alternative_fates",
    "FM07": "NAD_redox_balance",
}

BIOLOGICAL_INTERPRETATION = {
    "FM01": "Highly conserved central glycolytic framework",
    "FM02": "Conserved in Lactiplantibacillus but comparatively variable among LAB lineages",
    "FM03": "Highly conserved phosphoketolase/pentose-phosphate-associated framework",
    "FM04": "Most conserved module across sampled LAB lineages",
    "FM05": "Generally conserved with moderate lineage-associated variation",
    "FM06": "Major lineage-associated variation in alternative pyruvate fates",
    "FM07": "Broadly conserved redox-balancing framework with moderate variation",
}


def split_ids(value, prefix):
    """Extract identifiers such as GC001 or FM01 from catalogue cells."""
    if pd.isna(value):
        return []
    pattern = rf"{prefix}\d+"
    return re.findall(pattern, str(value))


def marker_status(prevalence):
    """Dataset-specific operational prevalence categories."""
    if prevalence >= CORE_THRESHOLD:
        return "CORE"
    if prevalence >= NEAR_CORE_THRESHOLD:
        return "NEAR_CORE"
    if prevalence >= VARIABLE_THRESHOLD:
        return "VARIABLE"
    return "RARE"


def save(df, filename):
    path = OUTDIR / filename
    df.to_csv(path, index=False)
    print(f"WROTE: {path.relative_to(ROOT)}  shape={df.shape}")


# ------------------------------------------------------------
# Load input data
# ------------------------------------------------------------

meta = pd.read_csv(METADATA)
unique = pd.read_csv(UNIQUE_MARKERS)
gm = pd.read_csv(GENE_MATRIX)

# ------------------------------------------------------------
# Basic integrity checks
# ------------------------------------------------------------

if meta["Genome_ID"].duplicated().any():
    raise ValueError("Duplicate Genome_ID values found in genome_metadata.csv")

if gm["Genome_ID"].duplicated().any():
    raise ValueError("Duplicate Genome_ID values found in fermentation_gene_matrix.csv")

meta_ids = set(meta["Genome_ID"])
gm_ids = set(gm["Genome_ID"])

if meta_ids != gm_ids:
    missing_in_matrix = sorted(meta_ids - gm_ids)
    missing_in_meta = sorted(gm_ids - meta_ids)
    raise ValueError(
        "Genome_ID mismatch between metadata and gene matrix.\n"
        f"Missing in matrix: {missing_in_matrix[:10]}\n"
        f"Missing in metadata: {missing_in_meta[:10]}"
    )

if len(meta) != 507 or len(gm) != 507:
    raise ValueError(
        f"Expected 507 genomes; metadata={len(meta)}, gene_matrix={len(gm)}"
    )

# ------------------------------------------------------------
# Collapse 31 catalogue assignments into 28 unique markers
# ------------------------------------------------------------

presence = gm[["Genome_ID"]].copy()

marker_gene_ids = {}
marker_modules = {}
marker_symbols = {}

for _, row in unique.iterrows():
    marker = row["Marker_ID"]
    gene_ids = split_ids(row["Gene_IDs"], "GC")
    module_ids = split_ids(row["Module_IDs"], "FM")

    if not gene_ids:
        raise ValueError(f"No Gene_ID found for {marker}")

    missing_cols = [g for g in gene_ids if g not in gm.columns]
    if missing_cols:
        raise ValueError(
            f"{marker}: gene-matrix columns missing: {missing_cols}"
        )

    presence[marker] = gm[gene_ids].max(axis=1).astype(int)

    marker_gene_ids[marker] = gene_ids
    marker_modules[marker] = module_ids
    marker_symbols[marker] = row["Gene_symbol"]

if len(marker_gene_ids) != 28:
    raise ValueError(
        f"Expected 28 unique biological markers; found {len(marker_gene_ids)}"
    )

# Add authoritative metadata
analysis = presence.merge(
    meta[
        [
            "Genome_ID",
            "Species_current",
            "Genus",
            "Accession",
            "Study_group",
        ]
    ],
    on="Genome_ID",
    how="left",
    validate="one_to_one",
)

marker_cols = list(marker_gene_ids.keys())

# Verify manuscript collection counts
group_counts = analysis["Study_group"].value_counts().to_dict()
expected_groups = {
    "core_LAB": 486,
    "fermentative_non_LAB": 10,
    "external_comparator": 8,
    "Enterococcus": 3,
}

for group, expected in expected_groups.items():
    observed = group_counts.get(group, 0)
    if observed != expected:
        raise ValueError(
            f"Unexpected count for {group}: observed={observed}, expected={expected}"
        )

core_lab = analysis[analysis["Study_group"] == "core_LAB"].copy()
lpb = core_lab[core_lab["Genus"] == "Lactiplantibacillus"].copy()

if len(lpb) != 462:
    raise ValueError(f"Expected 462 Lactiplantibacillus genomes; found {len(lpb)}")


# ============================================================
# 1. Lactiplantibacillus marker landscape
# ============================================================

rows = []

for marker in marker_cols:
    n_present = int(lpb[marker].sum())
    n_total = len(lpb)
    prevalence = n_present / n_total

    rows.append(
        {
            "Marker_ID": marker,
            "Gene_symbol": marker_symbols[marker],
            "Gene_IDs": ";".join(marker_gene_ids[marker]),
            "Module_IDs": ";".join(marker_modules[marker]),
            "N_present": n_present,
            "N_total": n_total,
            "Prevalence": prevalence,
            "Status": marker_status(prevalence),
        }
    )

lpb_landscape = pd.DataFrame(rows)
save(
    lpb_landscape,
    "Lactiplantibacillus_fermentation_landscape.csv",
)


# ============================================================
# 2. Lactiplantibacillus species-level marker landscape
# ============================================================

rows = []

species_counts = (
    lpb.groupby("Species_current")
    .size()
    .sort_values(ascending=False)
)

for species, n_genomes in species_counts.items():
    sub = lpb[lpb["Species_current"] == species]

    for marker in marker_cols:
        n_present = int(sub[marker].sum())
        prevalence = n_present / n_genomes

        rows.append(
            {
                "Species": species,
                "N_genomes": int(n_genomes),
                "Marker_ID": marker,
                "Gene_symbol": marker_symbols[marker],
                "N_present": n_present,
                "Prevalence": prevalence,
                "Status": marker_status(prevalence),
            }
        )

species_landscape = pd.DataFrame(rows)
save(
    species_landscape,
    "Lactiplantibacillus_species_marker_landscape.csv",
)


# ============================================================
# 3. Core-LAB genus × marker landscape
# ============================================================

rows = []

genus_counts = (
    core_lab.groupby("Genus")
    .size()
    .sort_values(ascending=False)
)

for genus, n_genomes in genus_counts.items():
    sub = core_lab[core_lab["Genus"] == genus]

    for marker in marker_cols:
        n_present = int(sub[marker].sum())
        prevalence = n_present / n_genomes

        rows.append(
            {
                "Genus": genus,
                "N_genomes": int(n_genomes),
                "Marker_ID": marker,
                "Gene_symbol": marker_symbols[marker],
                "N_present": n_present,
                "Prevalence": prevalence,
            }
        )

genus_marker = pd.DataFrame(rows)
save(
    genus_marker,
    "LAB_genus_fermentation_landscape.csv",
)


# ============================================================
# 4. Module membership using 28 unique biological markers
# ============================================================

module_markers = {m: [] for m in MODULE_NAMES}

for marker, modules in marker_modules.items():
    for module in modules:
        if module in module_markers:
            module_markers[module].append(marker)

expected_module_sizes = {
    "FM01": 10,
    "FM02": 2,
    "FM03": 6,
    "FM04": 2,
    "FM05": 2,
    "FM06": 4,
    "FM07": 5,
}

for module, expected in expected_module_sizes.items():
    observed = len(module_markers[module])
    if observed != expected:
        raise ValueError(
            f"{module}: expected {expected} unique markers; found {observed}"
        )


# ============================================================
# 5. Genus × module matrix
#
# Module value = mean prevalence of its unique biological
# markers within that genus.
# ============================================================

matrix_rows = []
long_rows = []

for genus, n_genomes in genus_counts.items():
    sub = core_lab[core_lab["Genus"] == genus]

    row = {
        "Genus": genus,
        "N_genomes": int(n_genomes),
    }

    for module in MODULE_NAMES:
        markers = module_markers[module]

        marker_prevalences = [
            sub[m].mean()
            for m in markers
        ]

        module_mean = sum(marker_prevalences) / len(marker_prevalences)
        row[module] = module_mean

        long_rows.append(
            {
                "Genus": genus,
                "N_genomes": int(n_genomes),
                "Module_ID": module,
                "Module_name": MODULE_NAMES[module],
                "N_unique_markers": len(markers),
                "Mean_marker_prevalence": module_mean,
            }
        )

    matrix_rows.append(row)

module_matrix = pd.DataFrame(matrix_rows)
module_matrix = module_matrix[
    ["Genus", "N_genomes"] + list(MODULE_NAMES.keys())
]

save(module_matrix, "LAB_genus_module_matrix.csv")

module_long = pd.DataFrame(long_rows)
save(module_long, "LAB_genus_module_completeness.csv")


# ============================================================
# 6. Inter-genus descriptive variability
#
# Genera are equally weighted.
# SD = sample standard deviation (ddof=1).
# Main manuscript comparison uses genera with N_genomes > 1.
# ============================================================

def variability_table(matrix_df):
    rows = []

    for module in MODULE_NAMES:
        values = matrix_df[module].astype(float)

        mean = values.mean()
        sd = values.std(ddof=1) if len(values) > 1 else 0.0
        value_range = values.max() - values.min()
        cv = sd / mean if mean != 0 else float("nan")

        lowest_idx = values.idxmin()
        highest_idx = values.idxmax()

        rows.append(
            {
                "Module_ID": module,
                "Module_name": MODULE_NAMES[module],
                "N_genera": len(values),
                "Mean_prevalence": mean,
                "SD": sd,
                "Range": value_range,
                "CV": cv,
                "Lowest_genus": matrix_df.loc[lowest_idx, "Genus"],
                "Highest_genus": matrix_df.loc[highest_idx, "Genus"],
            }
        )

    out = pd.DataFrame(rows)

    # Ranking defined in the project:
    # largest range first, SD used as secondary criterion.
    ranked = out.sort_values(
        ["Range", "SD"],
        ascending=[False, False],
    ).copy()

    rank_map = {
        module: rank
        for rank, module in enumerate(
            ranked["Module_ID"],
            start=1,
        )
    }

    out["Variability_rank"] = out["Module_ID"].map(rank_map)

    return out.sort_values("Variability_rank").reset_index(drop=True)


variability_all = variability_table(module_matrix)
save(
    variability_all,
    "LAB_module_intergenus_variability.csv",
)

multi_genus_matrix = module_matrix[module_matrix["N_genomes"] > 1].copy()

if len(multi_genus_matrix) != 4:
    raise ValueError(
        f"Expected 4 multi-genome core-LAB genera; found {len(multi_genus_matrix)}"
    )

variability_multi = variability_table(multi_genus_matrix)
save(
    variability_multi,
    "LAB_module_intergenus_variability_multigenome.csv",
)


# ============================================================
# 7. Marker prevalence by study group
# ============================================================

rows = []

for group, sub in analysis.groupby("Study_group", sort=True):
    n_genomes = len(sub)

    for marker in marker_cols:
        n_present = int(sub[marker].sum())

        rows.append(
            {
                "Study_group": group,
                "N_genomes": n_genomes,
                "Marker_ID": marker,
                "Gene_symbol": marker_symbols[marker],
                "N_present": n_present,
                "Prevalence": n_present / n_genomes,
            }
        )

group_marker = pd.DataFrame(rows)
save(
    group_marker,
    "fermentation_marker_group_prevalence.csv",
)


# ============================================================
# 8. Conservative LAB contrast
#
# core_LAB prevalence -
# max(external comparator prevalence,
#     fermentative non-LAB prevalence)
# ============================================================

rows = []

for marker in marker_cols:

    def prev(group):
        sub = analysis[analysis["Study_group"] == group]
        return float(sub[marker].mean())

    p_lab = prev("core_LAB")
    p_ext = prev("external_comparator")
    p_fnl = prev("fermentative_non_LAB")

    rows.append(
        {
            "Marker_ID": marker,
            "Gene_symbol": marker_symbols[marker],
            "core_LAB_prevalence": p_lab,
            "external_comparator_prevalence": p_ext,
            "fermentative_non_LAB_prevalence": p_fnl,
            "max_non_LAB_prevalence": max(p_ext, p_fnl),
            "LAB_contrast": p_lab - max(p_ext, p_fnl),
        }
    )

contrast = (
    pd.DataFrame(rows)
    .sort_values(
        ["LAB_contrast", "Marker_ID"],
        ascending=[False, True],
    )
    .reset_index(drop=True)
)

save(
    contrast,
    "fermentation_marker_LAB_contrast.csv",
)


# ============================================================
# 9. Module summary by study group
# ============================================================

rows = []

for group, sub in analysis.groupby("Study_group", sort=True):
    for module in MODULE_NAMES:
        markers = module_markers[module]

        mean_prevalence = sub[markers].mean().mean()

        genome_coverage = sub[markers].mean(axis=1)

        rows.append(
            {
                "Study_group": group,
                "N_genomes": len(sub),
                "Module_ID": module,
                "Module_name": MODULE_NAMES[module],
                "N_unique_markers": len(markers),
                "Mean_marker_prevalence": mean_prevalence,
                "Mean_genome_module_coverage": genome_coverage.mean(),
            }
        )

module_group = pd.DataFrame(rows)
save(
    module_group,
    "fermentation_module_group_summary.csv",
)


# ============================================================
# 10. Final manuscript module-results table
# ============================================================

lpb_status = lpb_landscape.set_index("Marker_ID")["Status"]

var_lookup = variability_multi.set_index("Module_ID")

rows = []

for module in MODULE_NAMES:
    markers = module_markers[module]

    statuses = [lpb_status[m] for m in markers]

    n_core = statuses.count("CORE")
    n_near = statuses.count("NEAR_CORE")

    if n_core == len(markers):
        module_status = "COMPLETE_CORE"
    elif n_core + n_near == len(markers):
        module_status = "NEAR_COMPLETE"
    else:
        module_status = "VARIABLE"

    v = var_lookup.loc[module]

    rows.append(
        {
            "Module_ID": module,
            "Module_name": MODULE_NAMES[module],
            "N_unique_markers": len(markers),
            "Lactiplantibacillus_CORE": n_core,
            "Lactiplantibacillus_NEAR_CORE": n_near,
            "Lactiplantibacillus_status": module_status,
            "Multi_genus_mean_prevalence": v["Mean_prevalence"],
            "Intergenus_SD": v["SD"],
            "Intergenus_range": v["Range"],
            "Intergenus_CV": v["CV"],
            "Variability_rank": int(v["Variability_rank"]),
            "Biological_interpretation": BIOLOGICAL_INTERPRETATION[module],
        }
    )

final_table = pd.DataFrame(rows)
save(
    final_table,
    "final_fermentation_module_results.csv",
)


# ============================================================
# Final audit summary
# ============================================================

print()
print("=" * 68)
print("REPRODUCIBILITY AUDIT")
print("=" * 68)

print(f"Total genomes:                  {len(analysis)}")
print(f"Core LAB genomes:               {len(core_lab)}")
print(f"Lactiplantibacillus genomes:    {len(lpb)}")
print(f"Unique biological markers:      {len(marker_cols)}")
print(f"Multi-genome core-LAB genera:   {len(multi_genus_matrix)}")

core_count = int((lpb_landscape["Status"] == "CORE").sum())
near_count = int((lpb_landscape["Status"] == "NEAR_CORE").sum())

print(f"Lactiplantibacillus CORE:       {core_count}/28")
print(f"Lactiplantibacillus NEAR_CORE:  {near_count}/28")

pfl = lpb_landscape[
    lpb_landscape["Gene_symbol"] == "pfl"
].iloc[0]

print(
    "pfl prevalence:                 "
    f"{int(pfl['N_present'])}/{int(pfl['N_total'])} "
    f"({pfl['Prevalence']:.4f})"
)

print()
print("Multi-genome genus module means:")
for module in MODULE_NAMES:
    mean_value = multi_genus_matrix[module].mean()
    print(f"  {module}: {mean_value:.3f}")

print()
print("Inter-genus variability ranking:")
for _, row in variability_multi.sort_values("Variability_rank").iterrows():
    print(
        f"  {int(row['Variability_rank'])}. {row['Module_ID']} "
        f"mean={row['Mean_prevalence']:.3f} "
        f"SD={row['SD']:.4f} "
        f"range={row['Range']:.4f} "
        f"CV={row['CV']:.4f}"
    )

print()
print("DONE")
