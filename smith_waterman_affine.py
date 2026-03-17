from enum import IntEnum
import numpy as np
import blosum


class Trace(IntEnum):
    STOP = 0
    LEFT = 1
    UP = 2
    DIAGONAL = 3


def fasta_reader(sequence_file):
    lines = open(sequence_file).readlines()
    sequence_name_row = lines[0][1:]
    sequence = ""
    for line in lines[1:]:
        sequence += line.strip()
    return sequence_name_row.replace(" ", "").strip(), sequence.strip()


def smith_waterman_affine(seq1, seq2, score_matrix, gap_open, gap_extend):
    row = len(seq1) + 1
    col = len(seq2) + 1

    open_cost = gap_open + gap_extend

    # Generating the empty matrices for storing scores
    matrix_M = np.zeros(shape=(row, col), dtype=float)
    matrix_Ix = np.zeros(shape=(row, col), dtype=float)
    matrix_Iy = np.zeros(shape=(row, col), dtype=float)

    # Generating the empty matrices for tracing
    tracing_M = np.zeros(shape=(row, col), dtype=int)
    tracing_Ix = np.zeros(shape=(row, col), dtype=int)
    tracing_Iy = np.zeros(shape=(row, col), dtype=int)

    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)
    max_matrix = "M"

    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):

            # Calculating Ix scores (gap in seq2, move up)
            ix_from_M = matrix_M[i - 1, j] - open_cost
            ix_from_Ix = matrix_Ix[i - 1, j] - gap_extend
            ix_from_Iy = matrix_Iy[i - 1, j] - open_cost

            matrix_Ix[i, j] = max(0, ix_from_M, ix_from_Ix, ix_from_Iy)

            if matrix_Ix[i, j] == 0:
                tracing_Ix[i, j] = Trace.STOP
            elif matrix_Ix[i, j] == ix_from_Ix:
                tracing_Ix[i, j] = Trace.UP
            elif matrix_Ix[i, j] == ix_from_M:
                tracing_Ix[i, j] = Trace.DIAGONAL
            elif matrix_Ix[i, j] == ix_from_Iy:
                tracing_Ix[i, j] = Trace.LEFT

            # Calculating Iy scores (gap in seq1, move left)
            iy_from_M = matrix_M[i, j - 1] - open_cost
            iy_from_Iy = matrix_Iy[i, j - 1] - gap_extend
            iy_from_Ix = matrix_Ix[i, j - 1] - open_cost

            matrix_Iy[i, j] = max(0, iy_from_M, iy_from_Iy, iy_from_Ix)

            if matrix_Iy[i, j] == 0:
                tracing_Iy[i, j] = Trace.STOP
            elif matrix_Iy[i, j] == iy_from_Iy:
                tracing_Iy[i, j] = Trace.LEFT
            elif matrix_Iy[i, j] == iy_from_M:
                tracing_Iy[i, j] = Trace.DIAGONAL
            elif matrix_Iy[i, j] == iy_from_Ix:
                tracing_Iy[i, j] = Trace.UP

            # Calculating M scores (match/mismatch, diagonal)
            substitution = score_matrix[seq1[i - 1]][seq2[j - 1]]
            m_from_M = matrix_M[i - 1, j - 1] + substitution
            m_from_Ix = matrix_Ix[i - 1, j - 1] + substitution
            m_from_Iy = matrix_Iy[i - 1, j - 1] + substitution

            matrix_M[i, j] = max(0, m_from_M, m_from_Ix, m_from_Iy)

            if matrix_M[i, j] == 0:
                tracing_M[i, j] = Trace.STOP
            elif matrix_M[i, j] == m_from_M:
                tracing_M[i, j] = Trace.DIAGONAL
            elif matrix_M[i, j] == m_from_Ix:
                tracing_M[i, j] = Trace.UP
            elif matrix_M[i, j] == m_from_Iy:
                tracing_M[i, j] = Trace.LEFT

            # Tracking the cell with the maximum score
            for mat_name, val in [("M", matrix_M[i, j]),
                                   ("Ix", matrix_Ix[i, j]),
                                   ("Iy", matrix_Iy[i, j])]:
                if val >= max_score:
                    max_index = (i, j)
                    max_score = val
                    max_matrix = mat_name

    # Initialising the variables for tracing
    aligned_seq1 = ""
    aligned_seq2 = ""
    current_aligned_seq1 = ""
    current_aligned_seq2 = ""
    (max_i, max_j) = max_index
    current_matrix = max_matrix

    # Tracing and computing the pathway with the local alignment
    while True:
        if current_matrix == "M":
            if tracing_M[max_i, max_j] == Trace.STOP:
                break
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = seq2[max_j - 1]

            if tracing_M[max_i, max_j] == Trace.DIAGONAL:
                next_matrix = "M"
            elif tracing_M[max_i, max_j] == Trace.UP:
                next_matrix = "Ix"
            elif tracing_M[max_i, max_j] == Trace.LEFT:
                next_matrix = "Iy"

            max_i = max_i - 1
            max_j = max_j - 1
            current_matrix = next_matrix

        elif current_matrix == "Ix":
            if tracing_Ix[max_i, max_j] == Trace.STOP:
                break
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = '-'

            if tracing_Ix[max_i, max_j] == Trace.UP:
                next_matrix = "Ix"
            elif tracing_Ix[max_i, max_j] == Trace.DIAGONAL:
                next_matrix = "M"
            elif tracing_Ix[max_i, max_j] == Trace.LEFT:
                next_matrix = "Iy"

            max_i = max_i - 1
            current_matrix = next_matrix

        elif current_matrix == "Iy":
            if tracing_Iy[max_i, max_j] == Trace.STOP:
                break
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[max_j - 1]

            if tracing_Iy[max_i, max_j] == Trace.LEFT:
                next_matrix = "Iy"
            elif tracing_Iy[max_i, max_j] == Trace.DIAGONAL:
                next_matrix = "M"
            elif tracing_Iy[max_i, max_j] == Trace.UP:
                next_matrix = "Ix"

            max_j = max_j - 1
            current_matrix = next_matrix

        aligned_seq1 = aligned_seq1 + current_aligned_seq1
        aligned_seq2 = aligned_seq2 + current_aligned_seq2

    # Reversing the order of the sequences
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    # Computing alignment positions (1-based)
    end_i, end_j = max_index
    consumed1 = sum(1 for c in aligned_seq1 if c != '-')
    consumed2 = sum(1 for c in aligned_seq2 if c != '-')
    start_i = end_i - consumed1 + 1
    start_j = end_j - consumed2 + 1

    return max_score, aligned_seq1, aligned_seq2, (start_i, end_i), (start_j, end_j)


def print_blast_alignment(aligned_seq1, aligned_seq2, seq1_range, seq2_range, line_width=60):
    """Print alignment in BLAST-like format"""
    seq1_pos = seq1_range[0]
    seq2_pos = seq2_range[0]

    for block_start in range(0, len(aligned_seq1), line_width):
        block1 = aligned_seq1[block_start:block_start + line_width]
        block2 = aligned_seq2[block_start:block_start + line_width]

        # Building the midline
        midline = ""
        for a, b in zip(block1, block2):
            if a == '-' or b == '-':
                midline += " "
            elif a == b:
                midline += "|"
            elif score_matrix[a][b] > 0:
                midline += "+"
            else:
                midline += " "

        # Counting non-gap characters for position tracking
        seq1_consumed = sum(1 for c in block1 if c != '-')
        seq2_consumed = sum(1 for c in block2 if c != '-')

        seq1_end = seq1_pos + seq1_consumed - 1
        seq2_end = seq2_pos + seq2_consumed - 1

        print(f"Query  {seq1_pos:>6}  {block1}  {seq1_end}")
        print(f"               {midline}")
        print(f"Sbjct  {seq2_pos:>6}  {block2}  {seq2_end}")
        print()

        seq1_pos = seq1_end + 1
        seq2_pos = seq2_end + 1


if __name__ == "__main__":
    file_1_name, file_1 = fasta_reader("fasta/mouse_seq.fasta")
    file_2_name, file_2 = fasta_reader("fasta/thaliana_seq.fasta")

    score_matrix = blosum.BLOSUM(62)

    best_score, aln1, aln2, range1, range2 = smith_waterman_affine(
        file_1, file_2, score_matrix, gap_open=11, gap_extend=1
    )

    print(f"Best local alignment score: {best_score}")
    print(f"Seq1 region: {range1[0]}..{range1[1]}")
    print(f"Seq2 region: {range2[0]}..{range2[1]}")
    print(f"Alignment length: {len(aln1)}")
    print()
    print_blast_alignment(aln1, aln2, range1, range2)