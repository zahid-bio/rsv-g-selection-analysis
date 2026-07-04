FLU_REGIONS_LOCKED = "2026-07-04"

H3_ANTIGENIC_SITES = {
    "A": [122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144,
          145, 146, 150, 152, 168],
    "B": [128, 129, 155, 156, 157, 158, 159, 160, 163, 165, 186, 187, 188, 189,
          190, 192, 193, 194, 196, 197, 198],
    "C": [44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294,
          297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312],
    "D": [96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177,
          179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218,
          219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248],
    "E": [57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94,
          109, 260, 261, 262, 265],
}

H3_ANTIGENIC_ALL = sorted({r for residues in H3_ANTIGENIC_SITES.values() for r in residues})

H3_HA1_LENGTH = 329

RSV_REGIONS_LOCKED = "2026-07-04"
RSV_G_CCR = list(range(157, 199))
RSV_G_CX3C = list(range(182, 187))
RSV_G_LENGTH = 298

if __name__ == "__main__":
    print(f"H3 antigenic residues locked {FLU_REGIONS_LOCKED}: "
          f"{len(H3_ANTIGENIC_ALL)} residues across sites "
          f"{sorted(H3_ANTIGENIC_SITES)} within HA1 length {H3_HA1_LENGTH}.")
    frac = len(H3_ANTIGENIC_ALL) / H3_HA1_LENGTH
    print(f"Antigenic residues cover {len(H3_ANTIGENIC_ALL)}/{H3_HA1_LENGTH} "
          f"= {frac*100:.1f}% of HA1 (this is the by-chance baseline).")
