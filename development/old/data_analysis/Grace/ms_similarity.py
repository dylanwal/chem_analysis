
import chem_analysis as ca
from chem_analysis.utils.vector_similarity import *


def main():
    ms_TCB = ca.ms.MSSignal.from_file(r"TCB.csv")
    ms_TCB_sample = ca.ms.MSSignal.from_file(r"TCB_sample.csv")
    ms_C10 = ca.ms.MSSignal.from_file(r"C10.csv")

    ms_TCB, ms_TCB_sample, ms_C10 = ca.ms.unify_mz(ms_TCB, ms_TCB_sample, ms_C10)
    ms_TCB_sample.y[40] = 0

    scorers = [
        cosine_distance, earth_movers_distance, euclidean_distance, manhattan_distance, jaccard_similarity,
        pearson_correlation, hamming_distance, kl_divergence, dice_similarity,
        weighted_recall_score, dot_product_soft_norm, dtw_distance, frechet_distance,
        sliding_cosine_distance
    ]

    for scorer in scorers:
        score = scorer(ms_TCB.y, ms_TCB_sample.y)
        print(f"{scorer.__name__}: {score}")
    print()
    for scorer in scorers:
        score = scorer(ms_TCB.y, ms_C10.y)
        print(f"{scorer.__name__}: {score}")

    fig = ca.plotting.signal(ms_TCB, normalize=True)
    fig = ca.plotting.signal(ms_TCB_sample, fig=fig, normalize=True)
    fig = ca.plotting.signal(ms_C10, fig=fig, normalize=True)
    fig.show()


if __name__ == '__main__':
    main()
