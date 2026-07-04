import sys, csv, json, random
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import functional_regions as fr

N_PERM = 10000
SEED = 42


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    sites_csv, target = sys.argv[1], sys.argv[2]
    n_perm = int(sys.argv[3]) if len(sys.argv) > 3 else 10000

    if n_perm < 0:
        sys.exit("Error: Number of permutations cannot be negative.")

    if target == "flu":
        region = set(fr.H3_ANTIGENIC_ALL)
        domain_len = fr.H3_HA1_LENGTH
        locked = fr.FLU_REGIONS_LOCKED
        residue_col = "h3_ha1_number"
        region_name = "H3 antigenic sites A-E (union)"
    elif target in ("rsv_a", "rsv_b", "rsv_combined"):
        region = set(fr.RSV_G_CCR)
        domain_len = fr.RSV_G_LENGTH
        locked = fr.RSV_REGIONS_LOCKED
        residue_col = "a2_residue"
        region_name = "RSV G Central Conserved Region (A2 157-198)"
    elif target in ("rsv_a_cx3c", "rsv_b_cx3c", "rsv_combined_cx3c"):
        region = set(fr.RSV_G_CX3C)
        domain_len = fr.RSV_G_LENGTH
        locked = fr.RSV_REGIONS_LOCKED
        residue_col = "a2_residue"
        region_name = "RSV G CX3C motif (A2 182-186)"
    else:
        sys.exit(f"Unknown target: {target}")

    selected = []
    for r in csv.DictReader(open(sites_csv)):
        raw = r.get(residue_col, "")
        if raw not in ("", "None", None):
            selected.append(int(raw))
    selected = sorted(set(selected))
    n = len(selected)

    if n > domain_len:
        raise ValueError(f"Number of selected sites ({n}) cannot exceed domain length ({domain_len}).")

    if n == 0:
        print("Warning: No selected sites map into the tested domain. Reporting p-value = 1.0.")
        observed = 0.0
        hits = 0
        expected = 0.0
        p_value = 1.0
        enrichment = 1.0
        null = []
    else:
        hits = sum(1 for s in selected if s in region)
        observed = hits / n

        random.seed(SEED)
        domain = list(range(1, domain_len + 1))
        null = []
        for _ in range(n_perm):
            draw = random.sample(domain, n)
            null.append(sum(1 for s in draw if s in region) / n)

        expected = sum(null) / n_perm if n_perm > 0 else 0.0
        ge = sum(1 for f in null if f >= observed)
        p_value = (ge + 1) / (n_perm + 1) if n_perm > 0 else 1.0
        enrichment = observed / expected if expected else float("inf")

    out = {
        "target": target,
        "region_name": region_name,
        "regions_locked": locked,
        "n_selected_in_domain": n,
        "selected_sites": selected,
        "hits_in_region": hits,
        "observed_fraction": observed,
        "expected_fraction": expected,
        "enrichment_ratio": enrichment,
        "p_value": p_value,
        "n_permutations": n_perm,
        "seed": SEED,
        "region_size": len(region),
        "domain_length": domain_len,
        "null_distribution": null,
    }
    outpath = Path("results") / f"{target}_enrichment.json"
    json.dump(out, open(outpath, "w"))

    print(f"target                    : {target} (regions locked {locked})")
    print(f"region                    : {region_name}")
    print(f"selected sites in domain  : {n}  -> {selected}")
    print(f"in region                 : {hits}")
    print(f"OBSERVED fraction         : {observed*100:.1f}%")
    print(f"EXPECTED by chance        : {expected*100:.1f}%  (region is "
          f"{len(region)}/{domain_len} = {len(region)/domain_len*100:.1f}% of domain)")
    print(f"ENRICHMENT ratio          : {enrichment:.2f}x")
    print(f"P-VALUE (one-sided)       : {p_value:.4f}")
    print("RESULT:", "SIGNIFICANT enrichment" if p_value < 0.05
          else "NOT significant")
    print(f"wrote -> {outpath}")


if __name__ == "__main__":
    main()
