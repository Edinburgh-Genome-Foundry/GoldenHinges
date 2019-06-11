from goldenhinges import OverhangsSelector

possible_overhangs_sets = [
    ["CCCT", "AACG", "ATCG", "GCTG", "TACA", "GAGT", "CCGA"],
    ["AACC", "TACA", "TAGA", "ATGC", "GATA", "CTCC", "GTAA"],
    ["AGTG", "CAGG", "ACTC", "AAAA", "AGAC", "CGAA", "ATAG"],
]


for possible_overhangs in possible_overhangs_sets:
    selector = OverhangsSelector(
        possible_overhangs=possible_overhangs, differences=2
    )
    compatible_overhangs = selector.generate_overhangs_set()
    if len(compatible_overhangs) > 6:

        print(
            "The following set has more than 6 overhangs:",
            compatible_overhangs,
        )
        break
