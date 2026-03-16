def read_fasta_one_sequence(path: str) -> str:
    seq_lines = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # If we already collected a sequence, stop at next header
                if seq_lines:
                    break
                continue
            # Accept sequence lines
            seq_lines.append(line)

    if not seq_lines:
        raise ValueError(f"No sequence found in {path}")

    return "".join(seq_lines).replace(" ", "").upper()


def read_scoring_matrix(path: str):
    # Read non-empty, non-comment lines
    lines = []
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            lines.append(ln)

    if not lines:
        raise ValueError(f"Matrix file '{path}' is empty or only comments.")

    # Header: first line should be the column symbols
    header = lines[0].split()
    cols = header

    score = {a: {} for a in cols}

    # Rows: each subsequent line: RowSymbol then scores for each column
    for ln in lines[1:]:
        toks = ln.split()
        if len(toks) < len(cols) + 1:
            continue
        row = toks[0]
        vals = toks[1:len(cols) + 1]
        if row not in score:
            score[row] = {}
        for c, v in zip(cols, vals):
            score[row][c] = int(v)

    # Make symmetric fallback (BLOSUM62 is symmetric)
    for a in list(score.keys()):
        for b, v in score[a].items():
            score.setdefault(b, {})
            score[b].setdefault(a, v)

    return score

def sw_affine_local_alignment(s1: str, s2: str, score, gap_open: int, gap_extend: int):
    n, m = len(s1), len(s2)

    open_cost = gap_open + gap_extend  # first gap char cost
    ext_cost  = gap_extend            # additional gap chars

    # DP matrices
    M  = [[0] * (m + 1) for _ in range(n + 1)]
    Ix = [[0] * (m + 1) for _ in range(n + 1)]
    Iy = [[0] * (m + 1) for _ in range(n + 1)]

    # Pointers: store (prev_matrix, prev_i, prev_j) or None if reset to 0
    PM  = [[None] * (m + 1) for _ in range(n + 1)]
    PIx = [[None] * (m + 1) for _ in range(n + 1)]
    PIy = [[None] * (m + 1) for _ in range(n + 1)]

    best_score = 0
    best_pos = (0, 0, "M")  # (i, j, matrix)

    for i in range(1, n + 1):
        a = s1[i - 1]
        for j in range(1, m + 1):
            b = s2[j - 1]

            # Ix: gap in s2 (move up: i-1, j)
            candidates_ix = [
                (0, None),
                (M[i - 1][j]  - open_cost, ("M",  i - 1, j)),
                (Ix[i - 1][j] - ext_cost,  ("Ix", i - 1, j)),
                (Iy[i - 1][j] - open_cost, ("Iy", i - 1, j)),
            ]
            Ix[i][j], PIx[i][j] = max(candidates_ix, key=lambda x: x[0])

            # Iy: gap in s1 (move left: i, j-1)
            candidates_iy = [
                (0, None),
                (M[i][j - 1]  - open_cost, ("M",  i, j - 1)),
                (Iy[i][j - 1] - ext_cost,  ("Iy", i, j - 1)),
                (Ix[i][j - 1] - open_cost, ("Ix", i, j - 1)),
            ]
            Iy[i][j], PIy[i][j] = max(candidates_iy, key=lambda x: x[0])

            # M: match/mismatch (diag: i-1, j-1)
            sub = score[a][b]
            candidates_m = [
                (0, None),
                (M[i - 1][j - 1]  + sub, ("M",  i - 1, j - 1)),
                (Ix[i - 1][j - 1] + sub, ("Ix", i - 1, j - 1)),
                (Iy[i - 1][j - 1] + sub, ("Iy", i - 1, j - 1)),
            ]
            M[i][j], PM[i][j] = max(candidates_m, key=lambda x: x[0])

            # track best among all three
            for mat, val in (("M", M[i][j]), ("Ix", Ix[i][j]), ("Iy", Iy[i][j])):
                if val > best_score:
                    best_score = val
                    best_pos = (i, j, mat)

    # Traceback from best_pos until score 0
    i, j, mat = best_pos
    aln1, aln2 = [], []

    def get_score(ii, jj, mm):
        return M[ii][jj] if mm == "M" else (Ix[ii][jj] if mm == "Ix" else Iy[ii][jj])

    while i > 0 and j > 0 and get_score(i, j, mat) > 0:
        if mat == "M":
            prev = PM[i][j]
            if prev is None:
                break
            aln1.append(s1[i - 1])
            aln2.append(s2[j - 1])
            mat, i, j = prev
        elif mat == "Ix":
            prev = PIx[i][j]
            if prev is None:
                break
            aln1.append(s1[i - 1])
            aln2.append("-")
            mat, i, j = prev
        else:  # Iy
            prev = PIy[i][j]
            if prev is None:
                break
            aln1.append("-")
            aln2.append(s2[j - 1])
            mat, i, j = prev

    aln1.reverse()
    aln2.reverse()
    aligned_s1 = "".join(aln1)
    aligned_s2 = "".join(aln2)

    # Compute start/end (1-based inclusive) using best end indices and consumed residues
    end_i, end_j, _ = best_pos
    consumed1 = sum(1 for c in aligned_s1 if c != "-")
    consumed2 = sum(1 for c in aligned_s2 if c != "-")
    start_i = end_i - consumed1 + 1
    start_j = end_j - consumed2 + 1

    return best_score, aligned_s1, aligned_s2, (start_i, end_i), (start_j, end_j)


def print_blast_like(aln1: str, aln2: str, score, start1: int, start2: int, width: int = 60):
    pos1, pos2 = start1, start2

    for k in range(0, len(aln1), width):
        c1 = aln1[k:k+width]
        c2 = aln2[k:k+width]

        mid = []
        for a, b in zip(c1, c2):
            if a == "-" or b == "-":
                mid.append(" ")
            elif a == b:
                mid.append("|")
            else:
                mid.append("+" if score[a][b] > 0 else " ")
        mid = "".join(mid)

        c1_used = sum(1 for x in c1 if x != "-")
        c2_used = sum(1 for x in c2 if x != "-")
        end1 = pos1 + c1_used - 1 if c1_used else pos1 - 1
        end2 = pos2 + c2_used - 1 if c2_used else pos2 - 1

        print(f"Query {pos1:>6}  {c1}  {end1}")
        print(f"{'':12}  {mid}")
        print(f"Sbjct {pos2:>6}  {c2}  {end2}\n")

        pos1 += c1_used
        pos2 += c2_used


def main():
    gap_open = 11
    gap_extend = 1
    seq1_path = "fasta/mouse_seq.fasta"
    seq2_path = "fasta/thaliana_seq.fasta"
    matrix_path = "blossum62.txt"
    s1 = read_fasta_one_sequence(seq1_path)
    s2 = read_fasta_one_sequence(seq2_path)
    Score = read_scoring_matrix(matrix_path)

    best_score, aln1, aln2, (s1_start, s1_end), (s2_start, s2_end) = sw_affine_local_alignment(
        s1, s2, Score, gap_open=gap_open, gap_extend=gap_extend
    )

    print("Local alignment (Smith–Waterman, affine gaps)")
    print(f"Gap penalties: open={gap_open}, extend={gap_extend} (gap cost = open + extend*length)")
    print()
    print("Best local alignment score:", best_score)
    print(f"Seq1 region: {s1_start}..{s1_end}")
    print(f"Seq2 region: {s2_start}..{s2_end}")
    print(f"Alignment length: {len(aln1)}")
    print()
    print("Aligned seq1:")
    print(aln1)
    print("Aligned seq2:")
    print(aln2)
    print()
    print_blast_like(aln1, aln2, Score, s1_start, s2_start, width=60)


if __name__ == "__main__":
    main()
