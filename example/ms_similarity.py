
import chem_analysis as ca
from chem_analysis.utils.vector_similarity import *


def main():
    print("hi")
    # ms_TCB = ca.ms.MSSignal.from_file(r"./data/ms/TCB.csv")
    # ms_TCB_sample = ca.ms.MSSignal.from_file(r"./data/ms/TCB_sample.csv")
    #
    # ms_TCB, ms_TCB_sample = ca.ms.unify_mz(ms_TCB, ms_TCB_sample)

    # scorers = [
    #     cosine_distance, earth_movers_distance, euclidean_distance, manhattan_distance, jaccard_similarity,
    #     pearson_correlation, hamming_distance, mahalanobis_distance, kl_divergence, dice_similarity,
    #     weighted_recall_score, dot_product_soft_norm, dtw_distance
    # ]
    #
    # for scorer in scorers:
    #     score = scorer(ms_TCB.y, ms_TCB_sample.y)
    #     print(f"{scorer.__name__}: {score}")


if __name__ == '__main__':
    main()
