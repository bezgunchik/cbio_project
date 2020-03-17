# loc - index of the cell of interest
# values - an array of frequencies (split line from the matrix)
# This function calculates score for one cell in a matrix
def get_cell_score(loc, values):
    val = float(values[loc])
    sigma = 0
    for i in range(len(values)):
        if (i != loc):
            sigma += abs(values[i] - val)
    return sigma / (len(values) - 1)


# values - an array of frequencies (split line from the matrix)
# This function calculates scores for each cell and returns the difference between two highest scores in a row
def get_cpg_score(values):
    cell_scores = []
    for i in range(len(values)):
        cell_scores.append(get_cell_score(i, values))
    cell_scores.sort(reverse=True)
    return cell_scores[0] - cell_scores[1]

