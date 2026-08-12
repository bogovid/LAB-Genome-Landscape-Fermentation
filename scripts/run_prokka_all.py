#!/usr/bin/env python3

import csv
import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

PROJECT = Path(__file__).resolve().parents[2]

METADATA = PROJECT / "metadata" / "processed" / "genome_metadata.csv"
OUTROOT = PROJECT / "analysis" / "prokka_annotations"
LOGROOT = PROJECT / "logs" / "prokka"
SUMMARY = PROJECT / "results" / "tables" / "prokka_annotation_status.tsv"

JOBS = int(os.environ.get("PROKKA_JOBS", "2"))

OUTROOT.mkdir(parents=True, exist_ok=True)
LOGROOT.mkdir(parents=True, exist_ok=True)
SUMMARY.parent.mkdir(parents=True, exist_ok=True)


def annotate(row):
    genome_id = row["Genome_ID"]
    accession = row["Accession"]
    species = row["Species_current"]
    fasta = Path(row["Full_path_wsl"])

    outdir = OUTROOT / genome_id
    done_file = outdir / ".DONE"
    log_file = LOGROOT / f"{genome_id}.log"

    if done_file.exists():
        return genome_id, accession, species, "SKIPPED_DONE"

    if not fasta.is_file():
        return genome_id, accession, species, "FASTA_NOT_FOUND"

    outdir.mkdir(parents=True, exist_ok=True)

    parts = species.split()
    genus = parts[0] if parts else "Bacteria"
    species_name = "_".join(parts[1:]) if len(parts) > 1 else "sp"

    cmd = [
        "prokka",
        "--outdir", str(outdir),
        "--prefix", genome_id,
        "--locustag", genome_id.replace("_", ""),
        "--genus", genus,
        "--species", species_name,
        "--cpus", "1",
        "--force",
        str(fasta),
    ]

    with log_file.open("w", encoding="utf-8") as log:
        log.write("COMMAND:\n")
        log.write(" ".join(cmd) + "\n\n")
        log.flush()

        result = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if result.returncode == 0:
        done_file.touch()
        return genome_id, accession, species, "DONE"

    return genome_id, accession, species, f"FAILED_{result.returncode}"


def main():
    with METADATA.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    print(f"Ukupno genoma: {len(rows)}")
    print(f"Paralelnih Prokka poslova: {JOBS}")
    print()

    results = []

    with ThreadPoolExecutor(max_workers=JOBS) as executor:
        futures = [executor.submit(annotate, row) for row in rows]

        for completed, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)

            genome_id, accession, species, status = result

            print(
                f"[{completed:>3}/{len(rows)}] "
                f"{genome_id}  {accession}  {status}",
                flush=True,
            )

    results.sort()

    with SUMMARY.open("w", encoding="utf-8") as handle:
        handle.write("Genome_ID\tAccession\tSpecies\tStatus\n")

        for genome_id, accession, species, status in results:
            handle.write(
                f"{genome_id}\t{accession}\t{species}\t{status}\n"
            )

    print()
    print(f"Status tabela: {SUMMARY}")


if __name__ == "__main__":
    main()
