import numpy as np
from sklearn.metrics import adjusted_rand_score


def rand_index(labels_true, labels_pred):
    """
    Compute the Rand Index between two clustering results.
    
    Parameters:
        labels_true: Ground truth (correct) cluster labels.
        labels_pred: Cluster labels to evaluate.
    
    Returns:
        float: The Rand Index score.
    """
    return adjusted_rand_score(labels_true, labels_pred)


if __name__ == "__main__":
    # Example usage
    labels_true = [0, 0, 1, 1]
    labels_pred = [0, 0, 1, 1]
    score = rand_index(labels_true, labels_pred)
    print(f"Rand Index: {score}")