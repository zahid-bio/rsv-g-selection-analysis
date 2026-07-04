import csv, json, random
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
RES = PROJECT / "results"

MUCIN_START = 199
MUCIN_END = 298
DOMAIN = 298
N_PERM = 10000
SEED = 42


def main():
    rows = list(csv.DictReader(open(RES / "rsv_combined_selected_a2.csv")))
    sites = sorted({int(r["a2_residue"]) for r in rows if r["a2_residue"]})
    n = len(sites)
    mucin = set(range(MUCIN_START, MUCIN_END + 1))
    hits = sum(1 for s in sites if s in mucin)
    obs = hits / n

    random.seed(SEED)
    domain = list(range(1, DOMAIN + 1))
    null = []
    for _ in range(N_PERM):
        draw = random.sample(domain, n)
        null.append(sum(1 for s in draw if s in mucin) / n)
    exp = sum(null) / N_PERM
    ge = sum(1 for f in null if f >= obs)
    p_val = (ge + 1) / (N_PERM + 1)
    enrich = obs / exp if exp else float("inf")

    out = {
        "analysis": "POST-HOC exploratory — NOT pre-registered",
        "region_name": f"C-terminal mucin-like variable domain (A2 {MUCIN_START}-{MUCIN_END})",
        "region_size": len(mucin),
        "domain_length": DOMAIN,
        "n_sites_tested": n,
        "site_residues": sites,
        "hits_in_region": hits,
        "observed_fraction": obs,
        "expected_fraction": exp,
        "enrichment_ratio": enrich,
        "p_value": p_val,
        "n_permutations": N_PERM,
        "seed": SEED,
        "caveat": ("The boundary A2 199-298 was chosen after observing that all "
                   "selected sites fell C-terminal to CCR. This p-value is "
                   "descriptive, not a hypothesis test. Report only as an "
                   "EXPLORATORY finding and do not use it to support a primary "
                   "claim."),
    }
    outpath = RES / "rsv_combined_mucin_posthoc.json"
    json.dump(out, open(outpath, "w"), indent=2)

    print(f"POST-HOC region        : A2 {MUCIN_START}-{MUCIN_END}  "
          f"({len(mucin)}/{DOMAIN} = {len(mucin)/DOMAIN*100:.1f}% of domain)")
    print(f"sites tested           : {n} -> {sites}")
    print(f"observed in mucin      : {hits}/{n} = {obs*100:.1f}%")
    print(f"expected by chance     : {exp*100:.1f}%")
    print(f"enrichment ratio       : {enrich:.2f}x")
    print(f"p-value (one-sided)    : {p_val:.5f}")
    print(f"NOTE: post-hoc, label as EXPLORATORY in the paper.")
    print(f"wrote -> {outpath}")


if __name__ == "__main__":
    main()
